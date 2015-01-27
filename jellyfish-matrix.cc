#pragma GCC diagnostic ignored "-Wunused-local-typedefs"

#include <string>
#include <fstream>
#include <vector>
#include <limits>

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

  /* default arguments */
  int abundance_threshold = 1; // count >= this are considered present in a sample
  int low_prevalence_threshold = 0; // kmers where num samples with kmer >= this are kept
  int high_prevalence_threshold = std::numeric_limits<int>::max(); // kmers where num samples with kmer <= this are kept
  int num_threads = 1;

  /* manadatory arguments */
  std::vector<std::string> in_files;
  std::string out_file;

  po::options_description desc("Form a matrix from many .jf files (max is 64), which are required to be in the same order and have exactly the same k-mers. The matrix is binary (presence/absence) and is stored using the jellyfish hash table where the value is the binary concatenation of each component .jf file. You may specify a low presence where only k-mers which appear in a number of samples greater than or equal to this are kept and likewise for a high presence. Lastly an option is provided for what value should be considered presence.");
  desc.add_options()
    ("help,h", "Show help message")
    ("input,i", po::value<std::vector<std::string> >(&in_files)->multitoken()->required(), "Input .jf files (ordered)")
    ("output,o", po::value< std::string>(&out_file)->required(), "File to write to")
    ("abundance,a", po::value<int>(&abundance_threshold)->default_value(abundance_threshold), "Counts >= this are considered present in a sample")
    ("min", po::value<int>(&low_prevalence_threshold)->default_value(low_prevalence_threshold), "kmers which appear in num samples >= this are kept")
    ("max", po::value<int>(&high_prevalence_threshold)->default_value(high_prevalence_threshold), "kmers which appear in num samples <= this are kept")
    ("threads,t", po::value<int>(&num_threads)->default_value(num_threads), "num threads to use");

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

  int num_files = in_files.size();
  jellyfish::file_header master_header;

  /* ---------- Open jf files ----------- */
  std::vector<binary_reader *> readers;
  for(int i = 0; i < num_files; ++i) {
    std::cerr << "Opening " << in_files[i] << " ... ";
    std::ifstream *ifstream = new std::ifstream(in_files[i]); // scary pointers
    if(!ifstream->good()) {
      std::cerr << "Error opening " << in_files[i] << std::endl;
      exit(1);
    }
    jellyfish::file_header *header = new jellyfish::file_header; // scary pointers
    std::cerr << "reading header ... ";
    header->read(*ifstream);
    if(i == 0) {
      master_header = *header;
      mer_dna::k(master_header.key_len() / 2);
    }
    //std::cerr << "header says " << header->size() << "," << header->key_len() << " ";
    std::cerr << "making reader ... ";
    binary_reader *reader = new binary_reader(*ifstream, header);
    readers.push_back(reader);
    std::cerr << "done." << std::endl;
  }

  mer_hash hash(master_header.size(), master_header.key_len(), num_files, num_threads, master_header.max_reprobe());
  hash.ary()->matrix(master_header.matrix());

  /* ---------- create the dumper ----------------- */
  std::auto_ptr<jellyfish::dumper_t<mer_array> > dumper;
  // number of bytes needed to store n samples
  int out_counter_length = num_files % 8 == 0 ? num_files / 8 : ((int) num_files / 8) + 1;
  dumper.reset(new binary_dumper(out_counter_length, master_header.key_len(), num_threads, out_file.c_str(), &master_header));
  hash.dumper(dumper.get());

  /* ---------- Iterate k-mers ---------- */
  {
    boost::timer::auto_cpu_timer t(std::cerr, 2);
    std::cerr << "=== Iterating " << master_header.size() << " k-mers ===" << std::endl;
    for(uint64_t i = 0; i < master_header.size(); ++i) {
      uint64_t presence = 0;
      int present = 0;
      for(int n = 0; n < num_files; ++n) {
        readers[n]->next();
        int count = readers[n]->val(); 
        presence <<= 1; // shift left, idempotent for start of loop
        if(count >= abundance_threshold) {
          ++present; // increase number of present
          ++presence; // set LSB to 1
        }
      }
      if(present >= low_prevalence_threshold && present <= high_prevalence_threshold) {
        hash.add(readers[0]->key(), presence);    
      }
    }
  }

  // TODO update the header.size() to reflect the number of actual stored k-mers
  /* ---------- Dump to output ---------- */
  {
    boost::timer::auto_cpu_timer t(std::cerr, 2);
    std::cerr << "=== Dumping jellyfish hash to " << out_file << " ===" << std::endl;
    dumper->one_file(true);
    dumper->dump(hash.ary());
  }

  return 0;
}
