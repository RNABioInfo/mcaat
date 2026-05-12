#include "cas/cas_writer.h"
#include <chrono>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <filesystem>
#include <algorithm>
#include <cmath>

namespace fs = std::filesystem;

// ── helpers ──────────────────────────────────────────────────────────────────

static std::string timestamp_now() {
    auto now   = std::chrono::system_clock::now();
    auto tt    = std::chrono::system_clock::to_time_t(now);
    std::ostringstream ss;
    ss << std::put_time(std::localtime(&tt), "%Y-%m-%d %H:%M:%S");
    return ss.str();
}

static std::string direction_label(SearchDirection d) {
    return (d == SearchDirection::DOWNSTREAM) ? "downstream" : "upstream";
}

// Left-pad an integer field so columns line up nicely
static std::string lpad(int v, int width) {
    std::string s = std::to_string(v);
    if (static_cast<int>(s.size()) < width)
        s = std::string(width - s.size(), ' ') + s;
    return s;
}

// ── CasWriter ────────────────────────────────────────────────────────────────

CasWriter::CasWriter(const std::string& output_dir, int max_lines_per_file)
    : output_dir_(output_dir), max_lines_per_file_(max_lines_per_file),
      generated_ts_(timestamp_now()) {
    fs::create_directories(output_dir_);
    OpenNextFile();
}

CasWriter::~CasWriter() {
    CloseFile();
}

std::string CasWriter::CurrentFilePath() const {
    return (fs::path(output_dir_) /
            ("CAS_Systems_" + std::to_string(file_index_) + ".txt")).string();
}

void CasWriter::OpenNextFile() {
    CloseFile();
    out_.open(CurrentFilePath(), std::ios::out | std::ios::trunc);
    if (!out_.is_open()) {
        throw std::runtime_error("CasWriter: cannot open " + CurrentFilePath());
    }
    line_count_ = 0;
    WriteFileHeader();
}

void CasWriter::CloseFile() {
    if (out_.is_open()) out_.close();
}

void CasWriter::Emit(const std::string& line) {
    // Roll to a new file when we're close to the line limit.
    // We avoid splitting a system block mid-way by rolling at block boundaries
    // (WriteSystemBlock checks headroom before writing).
    out_ << line << '\n';
    ++line_count_;
}

void CasWriter::EmitBlank() {
    Emit("");
}

// ── File header ──────────────────────────────────────────────────────────────

void CasWriter::WriteFileHeader() {
    Emit("MCAAT  Cas Gene Detection  v2.0.0");
    Emit("Generated  " + generated_ts_);
    // Stats line is written as a placeholder; we don't know totals yet at construction.
    // It will be updated by a second pass only if needed — for now leave informative stub.
    EmitBlank();
}

// ── Public write entry ───────────────────────────────────────────────────────

int CasWriter::Write(
    const std::vector<std::pair<std::string, std::vector<CasCassette>>>& results)
{
    // First pass: count totals for the summary header
    for (const auto& [repeat_seq, cassettes] : results) {
        for (const auto& c : cassettes) {
            if (!c.genes.empty()) {
                ++total_arrays_;
                total_genes_ += static_cast<int>(c.genes.size());
            }
        }
    }

    // Re-open first file and write the full header now that we have totals
    CloseFile();
    file_index_ = 1;
    out_.open(CurrentFilePath(), std::ios::out | std::ios::trunc);
    line_count_ = 0;
    Emit("MCAAT  Cas Gene Detection  v2.0.0");
    Emit("Generated  " + generated_ts_);
    Emit("Systems    " + std::to_string(total_arrays_));
    Emit("Genes      " + std::to_string(total_genes_));
    EmitBlank();

    int array_number = 0;
    for (const auto& [repeat_seq, cassettes] : results) {
        ++array_number;
        for (const auto& cassette : cassettes) {
            if (cassette.genes.empty()) continue;
            WriteSystemBlock(repeat_seq, cassette, array_number);
        }
    }

    std::cout << "  ▸ Cas output: " << total_arrays_ << " systems, "
              << total_genes_ << " genes  ("
              << file_index_ << " file" << (file_index_ > 1 ? "s" : "") << ")\n";
    return systems_written_;
}

// ── System block ─────────────────────────────────────────────────────────────

void CasWriter::WriteSystemBlock(const std::string& repeat_seq,
                                  const CasCassette& cassette,
                                  int array_number)
{
    // Estimate lines this block will take so we can roll before writing
    int est_lines = 5;  // separator + header + repeat + blank + trailing blank
    for (const auto& gene : cassette.genes) {
        est_lines += 2;  // gene stats line + blank line
        int aa_len = static_cast<int>(gene.amino_acids.size());
        est_lines += (aa_len > 0) ? ((aa_len + AA_LINE_WIDTH - 1) / AA_LINE_WIDTH) : 0;
    }

    // If block won't fit and we've written something, roll to next file
    if (line_count_ > 5 && line_count_ + est_lines > max_lines_per_file_) {
        ++file_index_;
        OpenNextFile();
    }

    std::string sep(80, '=');
    Emit(sep);

    // Build header line with aligned fields
    std::ostringstream hdr;
    hdr << "Array  " << lpad(array_number, 4)
        << "   subtype  " << (cassette.detected_type.empty() ? "Unknown" : cassette.detected_type)
        << "   span  " << cassette.total_distance_bp << " bp"
        << "   direction  " << direction_label(cassette.direction);
    if (!cassette.stop_reason_code.empty() && cassette.stop_reason_code != "LIMIT_REACHED") {
        hdr << "   stopped  " << cassette.stop_reason_code;
    }
    Emit(hdr.str());

    // Repeat line
    std::string repeat_display = repeat_seq.empty() ? "(unknown)" : repeat_seq;
    Emit("Repeat     " + repeat_display);
    EmitBlank();

    // Genes
    for (const auto& gene : cassette.genes) {
        // Gene stats line
        std::ostringstream gs;
        int aa_len = static_cast<int>(gene.amino_acids.size());
        // Score display: 2 decimal places, e.g. "0.82"
        char score_buf[16];
        std::snprintf(score_buf, sizeof(score_buf), "%.2f", gene.normalized_score);

        gs << "  " << std::left << std::setw(12) << gene.gene_family
           << "score  " << score_buf
           << "   length  " << lpad(aa_len, 4) << " aa"
           << "   offset  " << lpad(gene.distance_from_repeat, 6) << " bp";
        if (gene.is_putative) gs << "   [putative]";
        Emit(gs.str());

        // Amino-acid sequence
        WriteAASequence(gene.amino_acids);
        EmitBlank();
    }

    ++systems_written_;
}

void CasWriter::WriteAASequence(const std::string& aa) {
    if (aa.empty()) {
        Emit("  (no sequence)");
        return;
    }
    for (size_t i = 0; i < aa.size(); i += AA_LINE_WIDTH) {
        Emit("  " + aa.substr(i, AA_LINE_WIDTH));
    }
}
