#ifndef SPMUTIL_H
#define SPMUTIL_H

#include"spmtrx/spmcrs.h"

namespace spm {
	class spmutil {
		public:
			static spmcrs make_crs(const char *tsvfile);
			static spmcrs npmi(spmcrs& a, double threshold);
	};
}
#endif /** SPMUTIL_H **/
