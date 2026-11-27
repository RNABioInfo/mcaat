#include "protospacer_detection/label_features.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace mcaat::proto {
namespace {

constexpr uint64_t missing = std::numeric_limits<uint64_t>::max();

bool dna(const std::string &sequence) {
    return !sequence.empty() && std::all_of(sequence.begin(), sequence.end(), [](char c) {
               return c == 'A' || c == 'C' || c == 'G' || c == 'T';
           });
}

// ── A reader for the label file ───────────────────────────────────────────────
//
// The file is written by the labelling pipeline, so the shapes it can hold are
// known. This covers the whole of JSON except exponents, which the file never
// carries, and keeps the module free of a dependency.
class Json {
  public:
    enum class Kind { Null, Bool, Number, String, Array, Object };

    Kind kind = Kind::Null;
    bool boolean = false;
    double number = 0;
    std::string text;
    std::vector<Json> items;
    std::vector<std::pair<std::string, Json>> fields;

    const Json *find(const std::string &name) const {
        for (const auto &f : fields)
            if (f.first == name)
                return &f.second;
        return nullptr;
    }
    std::string string_or(const std::string &name, const std::string &fallback = {}) const {
        const Json *j = find(name);
        return j && j->kind == Kind::String ? j->text : fallback;
    }
    size_t number_or(const std::string &name, size_t fallback = 0) const {
        const Json *j = find(name);
        return j && j->kind == Kind::Number ? size_t(j->number) : fallback;
    }

    static Json parse(const std::string &source) {
        size_t at = 0;
        Json value = parse_value(source, at);
        skip(source, at);
        return value;
    }

  private:
    static void skip(const std::string &s, size_t &at) {
        while (at < s.size() && (s[at] == ' ' || s[at] == '\t' || s[at] == '\n' || s[at] == '\r'))
            ++at;
    }
    [[noreturn]] static void fail(size_t at, const std::string &what) {
        throw std::runtime_error("label file: " + what + " at byte " + std::to_string(at));
    }
    static void expect(const std::string &s, size_t &at, char c) {
        skip(s, at);
        if (at >= s.size() || s[at] != c)
            fail(at, std::string("expected ") + c);
        ++at;
    }
    static std::string parse_string(const std::string &s, size_t &at) {
        expect(s, at, '"');
        std::string out;
        while (at < s.size() && s[at] != '"') {
            if (s[at] == '\\') {
                if (++at >= s.size())
                    fail(at, "truncated escape");
                switch (s[at]) {
                case 'n': out += '\n'; break;
                case 't': out += '\t'; break;
                case 'r': out += '\r'; break;
                case 'b': out += '\b'; break;
                case 'f': out += '\f'; break;
                case 'u': {
                    if (at + 4 >= s.size())
                        fail(at, "truncated \\u escape");
                    unsigned code = std::stoul(s.substr(at + 1, 4), nullptr, 16);
                    at += 4;
                    if (code < 0x80) {
                        out += char(code);
                    } else if (code < 0x800) {
                        out += char(0xC0 | (code >> 6));
                        out += char(0x80 | (code & 0x3F));
                    } else {
                        out += char(0xE0 | (code >> 12));
                        out += char(0x80 | ((code >> 6) & 0x3F));
                        out += char(0x80 | (code & 0x3F));
                    }
                    break;
                }
                default: out += s[at];
                }
            } else {
                out += s[at];
            }
            ++at;
        }
        if (at >= s.size())
            fail(at, "unterminated string");
        ++at;
        return out;
    }
    static Json parse_value(const std::string &s, size_t &at) {
        skip(s, at);
        if (at >= s.size())
            fail(at, "no value");
        Json value;
        if (s[at] == '{') {
            ++at;
            value.kind = Kind::Object;
            skip(s, at);
            if (at < s.size() && s[at] == '}') {
                ++at;
                return value;
            }
            while (true) {
                skip(s, at);
                std::string name = parse_string(s, at);
                expect(s, at, ':');
                value.fields.emplace_back(std::move(name), parse_value(s, at));
                skip(s, at);
                if (at < s.size() && s[at] == ',') {
                    ++at;
                    continue;
                }
                expect(s, at, '}');
                return value;
            }
        }
        if (s[at] == '[') {
            ++at;
            value.kind = Kind::Array;
            skip(s, at);
            if (at < s.size() && s[at] == ']') {
                ++at;
                return value;
            }
            while (true) {
                value.items.push_back(parse_value(s, at));
                skip(s, at);
                if (at < s.size() && s[at] == ',') {
                    ++at;
                    continue;
                }
                expect(s, at, ']');
                return value;
            }
        }
        if (s[at] == '"') {
            value.kind = Kind::String;
            value.text = parse_string(s, at);
            return value;
        }
        if (s.compare(at, 4, "true") == 0) {
            at += 4;
            value.kind = Kind::Bool;
            value.boolean = true;
            return value;
        }
        if (s.compare(at, 5, "false") == 0) {
            at += 5;
            value.kind = Kind::Bool;
            value.boolean = false;
            return value;
        }
        if (s.compare(at, 4, "null") == 0) {
            at += 4;
            value.kind = Kind::Null;
            return value;
        }
        size_t begin = at;
        if (at < s.size() && (s[at] == '-' || s[at] == '+'))
            ++at;
        while (at < s.size() && (std::isdigit(static_cast<unsigned char>(s[at])) || s[at] == '.' ||
                                 s[at] == 'e' || s[at] == 'E' || s[at] == '-' || s[at] == '+'))
            ++at;
        if (begin == at)
            fail(at, "unrecognised value");
        value.kind = Kind::Number;
        value.number = std::stod(s.substr(begin, at - begin));
        return value;
    }
};

std::string tsv(double value) {
    std::ostringstream out;
    if (value == std::floor(value) && std::fabs(value) < 1e15)
        out << std::setprecision(0) << std::fixed << value;
    else
        out << std::setprecision(6) << value;
    return out.str();
}

} // namespace

// ── Class 1 ──────────────────────────────────────────────────────────────────

const std::vector<LabelEntry> LabelReader::empty_{};

LabelReader::LabelReader(const std::filesystem::path &file) {
    std::ifstream in(file, std::ios::binary);
    if (!in)
        throw std::runtime_error("cannot open label file " + file.string());
    std::string source((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    Json root = Json::parse(source);
    if (root.kind != Json::Kind::Object)
        throw std::runtime_error("label file is not a JSON object: " + file.string());

    for (const auto &[name, body] : root.fields) {
        if (name == "_meta") {
            for (const auto &[field, value] : body.fields)
                if (value.kind == Json::Kind::String)
                    meta_[field] = value.text;
                else if (value.kind == Json::Kind::Number)
                    meta_[field] = tsv(value.number);
            continue;
        }
        if (!name.empty() && name.front() == '_')
            continue; // _index and anything else reserved
        if (body.kind != Json::Kind::Object)
            continue;
        auto &bucket = by_dataset_[name];
        for (const auto &[spacer_key, flanks] : body.fields) {
            if (flanks.kind != Json::Kind::Array)
                continue;
            for (const auto &item : flanks.items) {
                if (item.kind != Json::Kind::Object)
                    continue;
                LabelEntry entry;
                entry.spacer_key = spacer_key;
                entry.id = item.string_or("label");
                entry.left = item.string_or("left");
                entry.spacer = item.string_or("spacer");
                entry.right = item.string_or("right");
                entry.record = item.string_or("record");
                entry.strand = item.string_or("strand");
                if (const Json *region = item.find("region");
                    region && region->kind == Json::Kind::Array && region->items.size() == 2) {
                    entry.region_begin = size_t(region->items[0].number);
                    entry.region_end = size_t(region->items[1].number);
                }
                if (const Json *arrays = item.find("arrays");
                    arrays && arrays->kind == Json::Kind::Array)
                    for (const auto &a : arrays->items)
                        if (a.kind == Json::Kind::String)
                            entry.arrays.push_back(a.text);
                if (entry.spacer.empty())
                    continue;
                bucket.push_back(std::move(entry));
            }
        }
    }
}

bool LabelReader::has(const std::string &dataset) const {
    return by_dataset_.count(dataset) != 0;
}

const std::vector<LabelEntry> &LabelReader::entries(const std::string &dataset) const {
    auto it = by_dataset_.find(dataset);
    return it == by_dataset_.end() ? empty_ : it->second;
}

std::vector<std::string> LabelReader::datasets() const {
    std::vector<std::string> names;
    names.reserve(by_dataset_.size());
    for (const auto &[name, entries] : by_dataset_)
        names.push_back(name);
    return names;
}

size_t LabelReader::size() const {
    size_t n = 0;
    for (const auto &[name, entries] : by_dataset_)
        n += entries.size();
    return n;
}

// ── Class 2 ──────────────────────────────────────────────────────────────────

TopologyFeatures::TopologyFeatures(const SDBG &graph) : graph_(graph), k_(size_t(graph.k())) {}

TopologyFeatures::TopologyFeatures(const SDBG &graph, Options options)
    : graph_(graph), options_(options), k_(size_t(graph.k())) {}

std::string TopologyFeatures::label(uint64_t node) const {
    if (node >= graph_.size())
        return {};
    std::vector<uint8_t> encoded(k_);
    graph_.GetLabel(node, encoded.data());
    std::string s;
    for (auto c : encoded) {
        if (c < 1 || c > 4)
            return {};
        s += "ACGT"[c - 1];
    }
    return s;
}

uint64_t TopologyFeatures::find_node(const std::string &word) const {
    if (word.size() != k_ || !dna(word))
        return missing;
    std::vector<uint8_t> encoded;
    encoded.reserve(word.size());
    for (char c : word)
        encoded.push_back(uint8_t(std::string("ACGT").find(c) + 1));
    uint64_t n = graph_.IndexBinarySearch(encoded.data());
    return n != missing && label(n) == word ? n : missing;
}

std::vector<TopologyFeatures::Arc> TopologyFeatures::neighbors(uint64_t node, bool backwards) const {
    std::vector<Arc> arcs;
    if (node >= graph_.size())
        return arcs;
    auto add = [&](uint64_t edge, uint64_t to) {
        auto c = graph_.GetW(edge);
        if (c > 4)
            c -= 4;
        if (c >= 1 && c <= 4 && graph_.IsValidEdge(edge))
            arcs.push_back({to, edge, double(graph_.EdgeMultiplicity(edge))});
    };
    if (!backwards) {
        uint64_t e = node;
        do {
            auto c = graph_.GetW(e);
            if (c > 4)
                c -= 4;
            if (c >= 1 && c <= 4 && graph_.IsValidEdge(e))
                add(e, graph_.Forward(e));
            --e;
        } while (e != missing && !graph_.IsLastOrTip(e));
    } else {
        uint64_t first = graph_.Backward(node);
        if (first >= graph_.size())
            return arcs;
        auto c = graph_.GetW(first);
        unsigned ends = graph_.IsLastOrTip(first);
        add(first, uint64_t(graph_.GetLastIndex(first)));
        for (uint64_t e = first + 1; ends < 5 && e < graph_.size(); ++e) {
            ends += graph_.IsLastOrTip(e);
            auto w = graph_.GetW(e);
            if (w == c)
                break;
            if (w == c + 4)
                add(e, uint64_t(graph_.GetLastIndex(e)));
        }
    }
    std::sort(arcs.begin(), arcs.end(), [](const Arc &a, const Arc &b) { return a.node < b.node; });
    return arcs;
}

char TopologyFeatures::added_base(const Arc &arc, bool backwards) const {
    if (backwards) {
        auto s = label(arc.node);
        return s.empty() ? 'N' : s.front();
    }
    auto c = graph_.GetW(arc.edge);
    if (c > 4)
        c -= 4;
    return c >= 1 && c <= 4 ? "ACGT"[c - 1] : 'N';
}

double TopologyFeatures::median(std::vector<double> values) {
    if (values.empty())
        return 0.;
    const size_t half = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + half, values.end());
    double m = values[half];
    if (values.size() % 2 == 0)
        m = (m + *std::max_element(values.begin(), values.begin() + half)) / 2.;
    return m;
}

double TopologyFeatures::variation(const std::vector<double> &values) {
    if (values.size() < 2)
        return 0.;
    double sum = std::accumulate(values.begin(), values.end(), 0.);
    double mean = sum / double(values.size());
    if (mean <= 0.)
        return 0.;
    double spread = 0.;
    for (double v : values)
        spread += (v - mean) * (v - mean);
    return std::sqrt(spread / double(values.size() - 1)) / mean;
}

// What lies ahead of one neighbour: how far the unbranching run goes, whether
// it dies, the multiplicity along it, and the nodes reachable inside the
// lookahead, which exist only so that a bubble can be recognised.
TopologyFeatures::Probe TopologyFeatures::probe(uint64_t from, const Arc &arc, bool backwards) const {
    Probe result;
    std::vector<double> depths{arc.count};

    // the unbranching run, stopping at the next decision
    uint64_t node = arc.node;
    uint64_t previous = from;
    while (result.run < options_.lookahead) {
        auto onward = neighbors(node, backwards);
        if (onward.empty()) {
            result.tip = true;
            break;
        }
        if (onward.size() > 1)
            break;
        if (onward.front().node == previous)
            break;
        previous = node;
        node = onward.front().node;
        depths.push_back(onward.front().count);
        ++result.run;
    }

    // visited set, bounded, so the probe does not walk a node twice
    std::unordered_set<uint64_t> seen;
    std::vector<std::pair<uint64_t, size_t>> stack{{arc.node, 0}};
    seen.insert(arc.node);
    while (!stack.empty() && seen.size() < options_.probe_nodes) {
        auto [at, depth] = stack.back();
        stack.pop_back();
        if (depth >= options_.lookahead)
            continue;
        for (const auto &next : neighbors(at, backwards))
            if (seen.insert(next.node).second)
                stack.emplace_back(next.node, depth + 1);
    }
    result.minimum = depths.empty() ? 0. : *std::min_element(depths.begin(), depths.end());
    result.cv = variation(depths);
    return result;
}

void TopologyFeatures::walk(const std::string &dataset, const LabelEntry &entry,
                            const std::string &region, size_t start, bool backwards,
                            std::vector<BranchFeature> &out) {
    if (region.size() < k_)
        return;

    // Anchor on the first k-mer of the region that the graph holds, starting
    // at the spacer boundary and moving outward.
    //
    // The whole flank is scanned. The reference region is a different strain's
    // sequence, so stretches of it are absent from this sample's graph; where
    // the walk loses the path it keeps moving outward one base at a time and
    // picks the path up again at the next k-mer the graph holds. Divergence
    // therefore costs the labels inside the diverged stretch and nothing more.
    size_t position = start;
    uint64_t node = missing;
    std::vector<double> running;
    size_t walked = 0, acquired = 0;
    double origin = 0.;

    const auto in_range = [&](size_t at) { return at + k_ <= region.size(); };
    const auto advance = [&](size_t &at) {
        if (backwards) {
            if (at == 0)
                return false;
            --at;
            return true;
        }
        ++at;
        return in_range(at);
    };

    while (in_range(position)) {
        if (node == missing) {
            node = find_node(region.substr(position, k_));
            if (node == missing) {
                if (!advance(position))
                    break;
                continue;
            }
            if (acquired++ == 0) {
                ++stats_.anchors;
                ++stats_.sides_mapped;
                auto arcs = neighbors(node, backwards);
                for (const auto &a : arcs)
                    origin = std::max(origin, a.count);
            } else {
                ++stats_.reanchors;
            }
        }

        // the base the reference continues into
        bool have_next = backwards ? position > 0 : position + k_ < region.size();
        if (!have_next)
            break;
        const char wanted = backwards ? region[position - 1] : region[position + k_];

        auto arcs = neighbors(node, backwards);
        int chosen = -1;
        for (size_t i = 0; i < arcs.size(); ++i)
            if (added_base(arcs[i], backwards) == wanted) {
                chosen = int(i);
                break;
            }

        if (arcs.size() > 1) {
            ++stats_.decisions;
            if (chosen < 0)
                ++stats_.decisions_skipped;
            else {
                double total = 0.;
                for (const auto &a : arcs)
                    total += a.count;
                // The walk has no history at its first decision, so the
                // multiplicity at the spacer boundary stands in for the median.
                const double level = running.empty() ? origin : median(running);

                std::vector<Probe> probes(arcs.size());
                for (size_t i = 0; i < arcs.size(); ++i)
                    probes[i] = probe(node, arcs[i], backwards);

                for (size_t i = 0; i < arcs.size(); ++i) {
                    BranchFeature f;
                    f.dataset = dataset;
                    f.label = entry.id;
                    f.spacer_key = entry.spacer_key;
                    f.backwards = backwards;
                    f.walked = walked;
                    f.node = arcs[i].node;
                    f.edge = arcs[i].edge;
                    f.multiplicity = arcs[i].count;
                    f.target = int(i) == chosen ? 1 : 0;

                    f.multiplicity_share = total > 0. ? arcs[i].count / total : 0.;

                    const Probe &p = probes[i];
                    f.run_ahead = p.run;
                    f.run_ends_tip = p.tip ? 1 : 0;
                    f.run_depth_min_ratio = level > 0. ? p.minimum / level : 0.;
                    f.run_depth_cv = p.cv;
                    f.in_degree = neighbors(arcs[i].node, !backwards).size();

                    if (f.target)
                        ++stats_.positives;
                    else
                        ++stats_.negatives;
                    out.push_back(std::move(f));
                }
            }
        }

        if (chosen < 0) {
            // The reference leaves the graph here, which a mismatch against
            // this sample does. Drop the path and keep scanning outward; the
            // loop above acquires a node again at the next k-mer the graph
            // holds.
            node = missing;
            running.clear();
            if (!advance(position))
                break;
            continue;
        }

        running.push_back(arcs[size_t(chosen)].count);
        if (running.size() > 2 * k_)
            running.erase(running.begin());
        node = arcs[size_t(chosen)].node;
        if (!advance(position))
            break;
        ++walked;
        ++stats_.steps;
    }

    if (acquired == 0)
        ++stats_.sides_empty;
}

std::vector<BranchFeature> TopologyFeatures::collect(const std::string &dataset,
                                                     const LabelEntry &entry) {
    std::vector<BranchFeature> out;
    ++stats_.entries;
    const std::string region = entry.region();
    if (region.size() <= k_ || !dna(region)) {
        ++stats_.unmapped;
        return out;
    }

    const size_t before = stats_.sides_mapped;

    // Right of the spacer: forward from the node ending at the spacer's end.
    if (entry.spacer_end() >= k_)
        walk(dataset, entry, region, entry.spacer_end() - k_, false, out);
    // Left of the spacer: backward from the node beginning at the spacer's start.
    if (entry.spacer_begin() + k_ <= region.size())
        walk(dataset, entry, region, entry.spacer_begin(), true, out);

    const size_t sides = stats_.sides_mapped - before;
    if (sides == 0 || (options_.require_both_sides && sides < 2)) {
        ++stats_.unmapped;
        if (options_.require_both_sides && sides < 2)
            out.clear();
        return out;
    }
    ++stats_.mapped;
    return out;
}

// ── Class 3 ──────────────────────────────────────────────────────────────────

// Eight identity and diagnostic columns, then the six features. target is the
// last column before the features so that a reader of the table can see where
// the fitting set begins.
std::vector<std::string> FeatureWriter::columns() {
    return {"dataset", "label", "spacer_key", "side", "walked", "node", "edge",
            "multiplicity", "target",
            // the features
            "multiplicity_share", "run_ahead", "run_ends_tip",
            "run_depth_min_ratio", "run_depth_cv", "in_degree"};
}

FeatureWriter::FeatureWriter(const std::filesystem::path &file) : file_(file) {
    if (file_.has_parent_path() && !file_.parent_path().empty())
        std::filesystem::create_directories(file_.parent_path());
    out_.open(file_, std::ios::binary);
    if (!out_)
        throw std::runtime_error("cannot write " + file_.string());
    auto names = columns();
    for (size_t i = 0; i < names.size(); ++i)
        out_ << names[i] << (i + 1 == names.size() ? '\n' : '\t');
}

FeatureWriter::~FeatureWriter() {
    if (out_.is_open())
        out_.close();
}

void FeatureWriter::write(const BranchFeature &r) {
    out_ << r.dataset << '\t' << r.label << '\t' << r.spacer_key << '\t'
         << (r.backwards ? "left" : "right") << '\t' << r.walked << '\t' << r.node << '\t'
         << r.edge << '\t' << tsv(r.multiplicity) << '\t' << r.target << '\t'
         << tsv(r.multiplicity_share) << '\t' << r.run_ahead << '\t'
         << r.run_ends_tip << '\t' << tsv(r.run_depth_min_ratio) << '\t'
         << tsv(r.run_depth_cv) << '\t' << r.in_degree << '\n';
    ++rows_;
}

void FeatureWriter::write(const std::vector<BranchFeature> &rows) {
    for (const auto &r : rows)
        write(r);
}

void FeatureWriter::finish(const TopologyFeatures::Stats &stats, const std::string &graph,
                           const std::string &labels, const std::vector<std::string> &datasets) {
    out_.flush();
    out_.close();
    std::filesystem::path side = file_;
    side += ".summary.json";
    std::ofstream out(side, std::ios::binary);
    if (!out)
        return;
    out << "{\n \"table\": \"" << file_.filename().string() << "\",\n";
    out << " \"rows\": " << rows_ << ",\n";
    out << " \"graph\": \"" << graph << "\",\n";
    out << " \"labels\": \"" << labels << "\",\n";
    out << " \"datasets\": [";
    for (size_t i = 0; i < datasets.size(); ++i)
        out << (i ? ", " : "") << '"' << datasets[i] << '"';
    out << "],\n";
    out << " \"entries\": " << stats.entries << ",\n";
    out << " \"entries_mapped\": " << stats.mapped << ",\n";
    out << " \"entries_unmapped\": " << stats.unmapped << ",\n";
    out << " \"sides_mapped\": " << stats.sides_mapped << ",\n";
    out << " \"sides_empty\": " << stats.sides_empty << ",\n";
    out << " \"anchors\": " << stats.anchors << ",\n";
    out << " \"reanchors\": " << stats.reanchors << ",\n";
    out << " \"steps\": " << stats.steps << ",\n";
    out << " \"decisions\": " << stats.decisions << ",\n";
    out << " \"decisions_skipped\": " << stats.decisions_skipped << ",\n";
    out << " \"positives\": " << stats.positives << ",\n";
    out << " \"negatives\": " << stats.negatives << "\n}\n";
}

} // namespace mcaat::proto
