#pragma GCC diagnostic ignored "-Wunused-local-typedefs"

#include <string>
#include <fstream>
#include <vector>
#include <cmath>

#include <boost/timer/timer.hpp>
#include <boost/program_options.hpp>

#include <jellyfish/err.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/hash_counter.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/mer_overlap_sequence_parser.hpp>
#include <jellyfish/whole_sequence_parser.hpp>
#include <jellyfish/mer_iterator.hpp>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/jellyfish.hpp>
#include <jellyfish/large_hash_iterator.hpp>

using jellyfish::mer_dna;
namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  //bool canonical = false;
  std::string in_file;

  po::options_description desc("Count occurrence of haplotype profiles in a jellyfish matrix");
  desc.add_options()
    ("help,h", "Show help message")
    //("canonical,c", po::bool_switch(&canonical)->default_value(canonical), "Group 101 and 010, labeled as 010")
    ("input,i", po::value< std::string>(&in_file)->required(), ".jf matrix file");

  /* ---------- parse arguments ---------------- */
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if(vm.count("help")) {
      std::cout << desc;
      exit(0);
    } 
    po::notify(vm);
  } catch(po::error &e) {
    std::cerr << "ERROR:" << e.what() << std::endl << desc;
    exit(1);
  }

  /* ---------- Get header from jf file ----------- */
  std::ifstream in_file_stream(in_file);
  if(!in_file_stream.good()) {
    std::cerr << "Error opening " << in_file << std::endl;
    exit(1);
  }

  jellyfish::file_header header;
  header.read(in_file_stream);
  mer_dna::k(header.key_len() / 2);

  int matrix_width = header.val_len();
  uint64_t canonical_cutoff = 1 << (matrix_width - 1); // choose leading bit of 0 to be canonical representation
  uint64_t canonical_mask = (~((uint64_t) 0)) >> (64 - matrix_width); // all ones in rightmost matrix_width bits
  std::vector<uint64_t> counts(canonical_cutoff);

  {
    boost::timer::auto_cpu_timer t(std::cerr, 2);
    std::cerr << "=== Iterating " << header.size() << " kmers ===" << std::endl;
    binary_reader reader(in_file_stream, &header);
    while(reader.next()) {
      uint64_t val = reader.val();
      if(val >= canonical_cutoff) {
        val = (~val) & canonical_mask;
      }
      ++counts[val];
    }
  }

  std::cerr << "value\tfrequency" << std::endl;

  for(uint64_t i = 0; i < counts.size(); ++i) {
    if(counts[i] > 0) {
      std::cout << i << "\t" << counts[i] << std::endl;
    }
  }

  return 0;
}
