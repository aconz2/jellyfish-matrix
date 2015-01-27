#pragma GCC diagnostic ignored "-Wunused-local-typedefs"

#include <string>
#include <fstream>
#include <vector>
#include <set>

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
/* ========== main ========== */

int main(int argc, char *argv[]) {

  /* default values */
  int num_threads = 1;
  uint64_t hash_size = 0;

  /* manadatory arguments */
  std::string jf_file;

  po::options_description desc("take in a .jf matrix, line seperated reads on stdin, output for each read, the number of haplotypes it contains");
  desc.add_options()
    ("help,h", "Show help message")
    ("jf,j", po::value< std::string>(&jf_file)->required(), ".jf file")
    //("size,s", po::value<uint64_t>(&hash_size)->default_value(hash_size), "Size of hash")
    ("threads,t", po::value<int>(&num_threads)->default_value(num_threads), "Number of threads to use");

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
  std::ifstream jf_file_stream(jf_file);
  if(!jf_file_stream.good()) {
    std::cerr << "Error opening " << jf_file << std::endl;
    exit(1);
  }
  jellyfish::file_header header;
  header.read(jf_file_stream);

  std::cerr << header.size() << std::endl;

  unsigned int kmer_length = header.key_len() / 2;
  hash_size = hash_size == 0 ? header.size() : hash_size;
  mer_dna::k(kmer_length);
  mer_hash hash(hash_size, header.key_len(), header.val_len(), num_threads, header.max_reprobe());
  hash.ary()->matrix(header.matrix());

  /* ---------- load the hash table with the k-mers from given file ---------- */
  {
    boost::timer::auto_cpu_timer t(std::cerr, 2);
    std::cerr << "=== Loading hash ===" << std::endl;
    binary_reader reader(jf_file_stream, &header);
    while(reader.next()) {
      hash.add(reader.key(), reader.val());
    }
  }

  /* ---------- Count k-mers from input ---------- */
  {
    boost::timer::auto_cpu_timer t(std::cerr, 2);
    std::cerr << "=== Examine reads on stdin ===" << std::endl;
    std::cerr << "Unique haplotype\tnot found k-mer" << std::endl;
    mer_array* ary = hash.ary();
    uint64_t val;
    mer_dna mer;
    for(std::string line; std::getline(std::cin, line);) {
      std::set<uint64_t> vals;
      int not_found = 0;
      for(unsigned int i = 0; i < line.size() - kmer_length; ++i) {
        mer = line.substr(i, kmer_length);
        //std::string substr = line.substr(i, kmer_length);
        //mer = substr;
        if(ary->get_val_for_key(mer, &val)) {
          vals.emplace(val);
        } else {
          not_found++;
        }
      }
      std::cout << vals.size() << "\t" << not_found << std::endl;
    }
  }
  return 0;
}
