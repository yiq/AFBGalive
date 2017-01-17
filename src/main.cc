#include <config.h>
#include <iostream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "sqlite3.h"

int compress_output = 0;


int main(int argc, char **argv){
	std::string g1kTable   = "g1k";
	std::string exacTable  = "exac";

	int c;
	// If the -b argument not provided, assume build GRCh37
	std::string build      = "GRCh37";

	while((c = getopt(argc, argv, ":zb:")) != -1) {
		switch (c) {
			case 'z':
				compress_output = 1;
				break;
			case 'b':
				// If the build argument was provided, append it to the exac and 1kg table names
        		build = optarg;
       			break;
			case ':':
				break;
			case '?':
				break;
			default:
				break;
		}
	}

	// Append the build onto the base table names for g1k and exac
	g1kTable  = g1kTable  + "_" + build;
	exacTable = exacTable + "_" + build;

	sqlite3 *db;
	int rc = sqlite3_open("afdata/af.db", &db);
	if(rc != SQLITE_OK) {
		std::cerr<<"Unable to open database file af.db"<<std::endl;
		return 1;
	}


	const char *errmsg;

	sqlite3_stmt *g1k_af_stmt;
	std::string g1k_af_sql = "SELECT val FROM " + g1kTable + " WHERE chrom=? AND pos=? AND ref=? AND alt=?;";
	rc = sqlite3_prepare_v2(db, g1k_af_sql.c_str(), strlen(g1k_af_sql.c_str()), &g1k_af_stmt, &errmsg);
	if(rc != SQLITE_OK) {
		std::cerr<<"("<<rc<<") Unable to compile select for g1k statement: "<<errmsg<<std::endl;
		return 1;
	}

	sqlite3_stmt *exac_af_stmt;
	std::string exac_af_sql = "SELECT val FROM " + exacTable + " WHERE chrom=? AND pos=? AND ref=? AND alt=?;";
	rc = sqlite3_prepare_v2(db, exac_af_sql.c_str(), strlen(exac_af_sql.c_str()), &exac_af_stmt, &errmsg);
	if(rc != SQLITE_OK) {
		std::cerr<<"("<<rc<<") Unable to compile select for exac statement: "<<errmsg<<std::endl;
		return 1;
	}

	htsFile *in_fp = hts_open("-", "rz");
	htsFile *out_fp;
	if(compress_output) {
		out_fp = hts_open("-", "wz");
	}
	else {
		out_fp = hts_open("-", "w");
	}

	bcf_hdr_t *hdr = vcf_hdr_read(in_fp);
	bcf_hdr_t *out_hdr = bcf_hdr_dup(hdr);

	int ret;

	ret = bcf_hdr_append(out_hdr, "##INFO=<ID=BGAF_1KG,Number=A,Type=Float,Description=\"Background Allele Frequency in the 1000 Genomes Project\">");
	ret = bcf_hdr_append(out_hdr, "##INFO=<ID=BGAF_EXAC,Number=A,Type=Float,Description=\"Background Allele Frequency in the Exome Aggregation Consortium (ExAC) Project\">");
	
	bcf_hdr_write(out_fp, out_hdr);

	bcf1_t *rec = bcf_init1();

	while(vcf_read1(in_fp, hdr, rec) >= 0) {
		int32_t chrom = rec->rid;
		int32_t pos = rec->pos;

		bcf_unpack(rec, BCF_UN_ALL);

		int nafs = rec->n_allele-1;

		float *afs_1kg = (float *)malloc(sizeof(float) * nafs);
		float *afs_exac = (float *)malloc(sizeof(float) * nafs);


		for(size_t i=1; i<rec->n_allele; i++) {
			std::string key = std::to_string(chrom) + ":" + \
							  std::to_string(pos) + ":" + \
							  rec->d.allele[0] + ":" + \
							  rec->d.allele[i];
			
			//std::string AFStr1KG = rStore->stringValueForKeyInNamespace(key, "1KG");
			//std::string AFStrExac = rStore->stringValueForKeyInNamespace(key, "EXAC");
			
			std::string AFStr1KG = "";
			std::string AFStrExac = "";

			sqlite3_bind_int(g1k_af_stmt, 1, chrom);
			sqlite3_bind_int(g1k_af_stmt, 2, pos);
			sqlite3_bind_text(g1k_af_stmt, 3, rec->d.allele[0], strlen(rec->d.allele[0]), NULL);
			sqlite3_bind_text(g1k_af_stmt, 4, rec->d.allele[i], strlen(rec->d.allele[i]), NULL);
			rc = sqlite3_step(g1k_af_stmt);
			if (rc == SQLITE_ROW) {
				const unsigned char * val = sqlite3_column_text(g1k_af_stmt, 0);
				std::string sval((const char *)(val));
				AFStr1KG = sval;
			}
			
			sqlite3_bind_int(exac_af_stmt, 1, chrom);
			sqlite3_bind_int(exac_af_stmt, 2, pos);
			sqlite3_bind_text(exac_af_stmt, 3, rec->d.allele[0], strlen(rec->d.allele[0]), NULL);
			sqlite3_bind_text(exac_af_stmt, 4, rec->d.allele[i], strlen(rec->d.allele[i]), NULL);
			rc = sqlite3_step(exac_af_stmt);
			if (rc == SQLITE_ROW) {
				const unsigned char * val = sqlite3_column_text(exac_af_stmt, 0);
				std::string sval((const char *)(val));
				AFStrExac = sval;
			}

			afs_1kg[i-1] = AFStr1KG != "" ? std::stof(AFStr1KG) : 0;
			afs_exac[i-1] = AFStrExac != "" ? std::stof(AFStrExac) : 0;

			sqlite3_clear_bindings(g1k_af_stmt);
			sqlite3_clear_bindings(exac_af_stmt);
			sqlite3_reset(g1k_af_stmt);
			sqlite3_reset(exac_af_stmt);
		}

		ret = bcf_update_info_float(out_hdr, rec, "BGAF_1KG", afs_1kg, nafs);
		ret = bcf_update_info_float(out_hdr, rec, "BGAF_EXAC", afs_exac, nafs);

		vcf_write1(out_fp, out_hdr, rec);

		free(afs_1kg);
		free(afs_exac);

	}

	sqlite3_finalize(g1k_af_stmt);
	sqlite3_finalize(exac_af_stmt);

	bcf_hdr_destroy(out_hdr);
	bcf_hdr_destroy(hdr);
	hts_close(in_fp);
	hts_close(out_fp);

	sqlite3_close(db);
}
