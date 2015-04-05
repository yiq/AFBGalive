#include <iostream>
#include <string>
#include <cstdio>
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "RedisKVStore.h"

int main(int argc, const char *argv[]) {

	if(argc<2) {
		std::cerr<<"usage %s <namespace> < file"<<std::endl;
		return 1;
	}

	std::string ns(argv[1]);

	YiCppLib::RedisKVStore::pointer rStore(new YiCppLib::RedisKVStore("/opt/redis/var/run/redis.sock"));

	htsFile *fp = hts_open("-", "rg");
	bcf_hdr_t *hdr = vcf_hdr_read(fp);
	bcf1_t *rec = bcf_init1();

	int lastChrom = -1;

	while(vcf_read1(fp, hdr, rec) >= 0) {
		int32_t chrom = rec->rid;
		int32_t pos = rec->pos;

		if(lastChrom != chrom) {
			const char * name = bcf_hdr_id2name(hdr, chrom);
			std::cout<<std::endl;
			std::cout<<"chr" + std::string(name) + " "<<std::flush;
			lastChrom = chrom;
		}

		if(pos % 1000000 == 0) std::cout<<"#"<<std::flush;

		bcf_unpack(rec, BCF_UN_STR);
		//bcf_info_t *afInfo = bcf_get_info(hdr, rec, "AF");
		float *afs = new float[rec->n_allele-1];
		int nafs = rec->n_allele-1;

		int ret = bcf_get_info_float(hdr, rec, "AF", &afs, &nafs);

		for(size_t i=1; i<rec->n_allele; i++) {
			std::string key = std::to_string(chrom) + ":" + \
							  std::to_string(pos) + ":" + \
							  rec->d.allele[0] + ":" + \
							  rec->d.allele[i];

			std::string val = std::to_string(afs[i-1]);
			rStore->setStringValueForKeyInNamespace(val, key, ns);
		}
		delete [] afs;
	}
	std::cout<<std::endl;

	bcf_hdr_destroy(hdr);
	hts_close(fp);
}
