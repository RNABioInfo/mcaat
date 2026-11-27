// Collects the training table for a learned branch decision.
//
// For one sample, every reference-confirmed protospacer region in the label
// file is mapped into that sample's graph and walked outward from the spacer.
// At each node with more than one neighbour, one row is written per neighbour:
// the neighbour the reference continues into carries target=1, its siblings
// target=0. Every other column is read from topology and edge multiplicity, so
// nothing in the table depends on a reference being present.
//
//   label_features --graph PREFIX --labels labels_by_dataset.json
//                  --dataset SRR14017550 --output features.tsv
//                  [--lookahead 300] [--probe-nodes 2000] [--both-sides]
//                  [--list]
//
// Writes the table and features.tsv.summary.json beside it.
#include "protospacer_detection/label_features.h"
#include "sdbg/sdbg.h"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace mcaat::proto;

namespace {

void usage() {
    std::cerr << "label_features --graph PREFIX --labels FILE --dataset NAME --output FILE\n"
              << "               [--lookahead N] [--probe-nodes N] [--both-sides] [--list]\n";
}

} // namespace

int main(int argc, char **argv) {
    std::string graph_prefix, labels, dataset, output;
    TopologyFeatures::Options options;
    bool list = false;

    for (int i = 1; i < argc; ++i) {
        std::string flag = argv[i];
        auto value = [&]() -> std::string {
            if (i + 1 >= argc)
                throw std::invalid_argument("missing value for " + flag);
            return argv[++i];
        };
        if (flag == "--graph")
            graph_prefix = value();
        else if (flag == "--labels")
            labels = value();
        else if (flag == "--dataset")
            dataset = value();
        else if (flag == "--output")
            output = value();
        else if (flag == "--lookahead")
            options.lookahead = std::stoul(value());
        else if (flag == "--probe-nodes")
            options.probe_nodes = std::stoul(value());
        else if (flag == "--both-sides")
            options.require_both_sides = true;
        else if (flag == "--list")
            list = true;
        else if (flag == "-h" || flag == "--help") {
            usage();
            return 0;
        } else {
            std::cerr << "unknown flag " << flag << "\n";
            usage();
            return 2;
        }
    }

    if (labels.empty()) {
        usage();
        return 2;
    }

    try {
        LabelReader reader(labels);

        if (list) {
            std::cout << "datasets in " << labels << ": " << reader.datasets().size()
                      << ", entries: " << reader.size() << "\n";
            for (const auto &name : reader.datasets())
                std::cout << name << '\t' << reader.entries(name).size() << '\n';
            return 0;
        }

        if (graph_prefix.empty() || output.empty()) {
            usage();
            return 2;
        }
        if (dataset.empty()) {
            auto names = reader.datasets();
            if (names.size() != 1) {
                std::cerr << "--dataset is required: the label file holds " << names.size()
                          << " datasets and a graph belongs to one of them\n";
                return 2;
            }
            dataset = names.front();
        }
        if (!reader.has(dataset)) {
            std::cerr << "no labels for " << dataset << " in " << labels << "\n";
            return 1;
        }

        const auto &entries = reader.entries(dataset);
        std::cerr << "dataset " << dataset << ": " << entries.size() << " labelled regions\n";

        SDBG graph;
        graph.LoadFromFile(graph_prefix.c_str());
        std::cerr << "graph " << graph_prefix << ": k=" << graph.k() << ", " << graph.size()
                  << " edges\n";

        TopologyFeatures features(graph, options);
        FeatureWriter writer(output);

        const auto began = std::chrono::steady_clock::now();
        for (const auto &entry : entries)
            writer.write(features.collect(dataset, entry));
        const auto seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - began).count();

        const auto &s = features.stats();
        writer.finish(s, graph_prefix, labels, {dataset});

        std::cerr << "regions mapped " << s.mapped << " of " << s.entries << ", sides mapped "
                  << s.sides_mapped << ", empty " << s.sides_empty << "\n"
                  << "steps " << s.steps << ", re-anchors " << s.reanchors << "\n"
                  << "decisions " << s.decisions << ", of them without the reference base "
                  << s.decisions_skipped << "\n"
                  << "rows " << writer.rows() << " (positives " << s.positives << ", negatives "
                  << s.negatives << ")\n"
                  << "wrote " << output << " in " << seconds << " s\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "label_features: " << error.what() << "\n";
        return 1;
    }
}
