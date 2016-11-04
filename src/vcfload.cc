#include <iostream>
#include <string>
#include <cstdio>
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "sqlite3.h"

int main(int argc, const char *argv[]) {

	if(argc<3) {
		std::cerr<<"usage %s <db_filename> <tablename> < file"<<std::endl;
		return 1;
	}

	std::string ns(argv[2]);

	htsFile *fp = hts_open("-", "rg");
	bcf_hdr_t *hdr = vcf_hdr_read(fp);
	bcf1_t *rec = bcf_init1();

	sqlite3 *db;
	int rc = sqlite3_open(argv[1], &db);

	if(rc != SQLITE_OK) {
		std::cerr<<"Unable to open database file "<<argv[1]<<std::endl;
		return 1;
	}
	std::cerr<<"DB Connection Open"<<std::endl;

	std::string sql_stmt = "CREATE TABLE " + ns + " ( chrom, pos, ref, alt, val);";
	char *errmsg;
	rc = sqlite3_exec(db, sql_stmt.c_str(), NULL, NULL, &errmsg);

	if(rc != SQLITE_OK) {
		std::cerr<<"Unable to create database table "<<ns<<": "<<errmsg<<std::endl;
		return 1;
	}
	std::cerr<<"DB Table created"<<std::endl;

	sqlite3_exec(db, "PRAGMA synchronous = OFF", NULL, NULL, &errmsg);
	sqlite3_exec(db, "PRAGMA journal_mode = MEMORY", NULL, NULL, &errmsg);

	sqlite3_stmt *insertStmt;
	std::string insert_sql = "INSERT INTO " + ns + " (chrom, pos, ref, alt, val) values (?,?,?,?,?);";
	rc = sqlite3_prepare_v2(db, insert_sql.c_str(), strlen(insert_sql.c_str()), &insertStmt, NULL);
	if(rc != SQLITE_OK) {
		std::cerr<<"Unable to compile insert statement"<<std::endl;
		return 1;
	}

	std::cerr<<"DB stmt prepared"<<std::endl;


	int lastChrom = -1;
	int n = 0;

	sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errmsg);

	while(vcf_read1(fp, hdr, rec) >= 0) {
		int32_t chrom = rec->rid;
		int32_t pos = rec->pos;

		if(lastChrom != chrom) {
			sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &errmsg);
			const char * name = bcf_hdr_id2name(hdr, chrom);
			std::cout<<std::endl;
			std::cout<<"chr" + std::string(name) + " "<<std::flush;
			lastChrom = chrom;
			sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errmsg);
		}

		if(n > 0 && n % 50000 == 0) std::cout<<"#"<<std::flush;

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

			//std::cout << "val=" << val.c_str() << std::endl;

			sqlite3_bind_int(insertStmt, 1, chrom);
			sqlite3_bind_int(insertStmt, 2, pos);
			sqlite3_bind_text(insertStmt, 3, rec->d.allele[0], strlen(rec->d.allele[0]), 0);
			sqlite3_bind_text(insertStmt, 4, rec->d.allele[i], strlen(rec->d.allele[i]), 0);
			sqlite3_bind_text(insertStmt, 5, val.c_str(), strlen(val.c_str()), 0);
			sqlite3_step(insertStmt);
			sqlite3_clear_bindings(insertStmt);
			sqlite3_reset(insertStmt);
			n++;

		}
		delete [] afs;
	}

	sqlite3_finalize(insertStmt);
	std::cout<<std::endl;

	sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &errmsg);


	bcf_hdr_destroy(hdr);
	hts_close(fp);
	sqlite3_close(db);
}
