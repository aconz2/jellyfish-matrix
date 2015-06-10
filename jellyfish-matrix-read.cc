#pragma GCC diagnostic ignored "-Wunused-local-typedefs"

#include <string>
#include <fstream>
#include <iterator>
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

struct fastq {
  std::string id;
  std::string seq;
  std::string plus;
  std::string qual;

  friend std::istream & operator>>(std::istream &is, fastq &record) {
    std::getline(is, record.id);
    std::getline(is, record.seq);
    std::getline(is, record.plus);
    std::getline(is, record.qual);
    return is;
  };
};

struct fasta {
  std::string id;
  std::string seq;

  friend std::istream & operator>>(std::istream &is, fasta &record) {
    std::getline(is, record.id);
    std::getline(is, record.seq);
    return is;
  };

};

using jellyfish::mer_dna;
namespace po = boost::program_options;
/* ========== main ========== */

int main(int argc, char *argv[]) {

  /* default values */
  int num_threads = 1;

  /* manadatory arguments */
  std::string jf_file;
  std::string fasta_file;
  std::string fastq_file;

  po::options_description desc("take in a .jf matrix, fastq file on stdin, output for each read, the number of haplotypes it contains");
  desc.add_options()
    ("help,h", "Show help message")
    ("jf,j", po::value<std::string>(&jf_file), ".jf file")
    ("fasta", po::value<std::string>(&fasta_file), "fasta file")
    ("fastq", po::value<std::string>(&fastq_file), "fastq file");
    //("threads,t", po::value<int>(&num_threads)->default_value(num_threads), "Number of threads to use")

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

  if(fasta_file == "" && fastq_file == "") {
    std::cerr << "ERROR: require one of fasta or fastq" << std::endl;
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

  uint64_t hash_size = header.size();
  unsigned int kmer_length = header.key_len() / 2;
  std::cerr << hash_size << std::endl;
  std::cerr << kmer_length << std::endl;

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

  /*{
    boost::timer::auto_cpu_timer t(std::cerr, 2);
    std::cerr << "=== Examine reads on stdin ===" << std::endl;
    mer_array* ary = hash.ary();
    uint64_t val;
    mer_dna mer;
    std::istream_iterator<fastq> fastq_iter(std::cin);
    std::istream_iterator<fastq> EOS;
    while(fastq_iter != EOS) {
      fastq record = *fastq_iter++;
      std::cout << record.id << std::endl;
      for(unsigned int i = 0; i < record.seq.size() - kmer_length + 1; ++i) {
        mer = record.seq.substr(i, kmer_length);
        if(!ary->get_val_for_key(mer, &val)) {
          std::cout << ".";
        } else {
          std::cout << val << " ";
        }
      }
      std::cout << std::endl;
    }
 }*/

  mer_array* ary = hash.ary();
  uint64_t val;
  mer_dna mer;

  // CODE DUPLICATION > 9000
  if(fasta_file != "") {
    std::ifstream fasta_if(fasta_file);
    std::istream_iterator<fasta> fasta_iter(fasta_if);
    std::istream_iterator<fasta> EOS;
    boost::timer::auto_cpu_timer t(std::cerr, 2);
    std::cerr << "=== Examine fasta from " << fasta_file <<  "===" << std::endl;
    while(fasta_iter != EOS) {
      fasta record = *fasta_iter++;
      std::cout << record.id << std::endl;
      for(unsigned int i = 0; i < record.seq.size() - kmer_length + 1; ++i) {
        // have to deal with 'N' values here (jellyfish deals with lowercases)
        if(!(mer.from_chars(record.seq.substr(i, kmer_length).c_str()))) {
          continue;
        }
        // and make sure we are getting the canonical
        if(ary->get_val_for_key(mer.get_canonical(), &val)) {
          std::cout << val << "\t";
        }
      }
      std::cout << std::endl;
    }
  } else {
    std::ifstream fastq_if(fastq_file);
    std::istream_iterator<fastq> fastq_iter(fastq_if);
    std::istream_iterator<fastq> EOS;
    boost::timer::auto_cpu_timer t(std::cerr, 2);
    std::cerr << "=== Examine fastq from " << fastq_file <<  "===" << std::endl;
    while(fastq_iter != EOS) {
      fastq record = *fastq_iter++;
      std::cout << record.id << std::endl;
      for(unsigned int i = 0; i < record.seq.size() - kmer_length + 1; ++i) {
        // have to deal with 'N' values here (jellyfish deals with lowercases)
        if(!(mer.from_chars(record.seq.substr(i, kmer_length).c_str()))) {
          continue;
        }
        // and make sure we are getting the canonical
        if(ary->get_val_for_key(mer.get_canonical(), &val)) {
          std::cout << val << "\t";
        }
      }
      std::cout << std::endl;
    }
  }
 
  return 0;
}
