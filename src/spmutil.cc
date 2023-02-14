#include"spmtrx/spmutil.h"
#include<string>
#include<fstream>
#include<algorithm>

using namespace spm;
using namespace std;

spmcrs spmutil::make_crs(const char *tsv) {
	ifstream f(tsv);
	if (!f) {
		throw "failed to open tsv";
	}
	string buf;
	vector<string> source;
	while (getline(f, buf)) {
		if (buf[buf.size()-1] == '\n') {
			buf[buf.size()-1] = '\0';
		}
		if (buf.find("\t") == string::npos) {
			continue;
		}
		source.push_back(buf);
	}
	spmcrs a;
	sort(source.begin(), source.end(), 
			[](const string& a, const string& b){
			return atoi(a.c_str()) < atoi(b.c_str());
			});
	int p_instance = 0;
	auto i = 0;
	for (auto it = source.begin(); it != source.end(); ++it) {
		int instance = 0;
		int pattern = 0;
		double score = 0.;
		sscanf(it->c_str(), "%d\t%d\t%lf", &instance, &pattern, &score);
		if (p_instance != instance) {
			a.cr();
			++i;
		}
		a.set(pattern, score);
		p_instance = instance;
	}
	a.cr();
	return a;
}

spmcrs spmutil::npmi(spmcrs& a, double threshold) {
	spmcrs npmi;
	spmcrs d = a.diag();
	double n = d.trace();
	spmcrs d_t = a.transpose().diag();

	auto i = 0;
	for (auto it = a.begin(); it != a.end(); ++it, ++i) {
		for (auto j = *it; j < *(it+1); ++j) {
			auto col = a.get_col(j);
			auto freq = a.get_val(j);
			double v = (log(freq)+log(n)-log(d.get_val(i, i))-log(d_t.get_val(col, col))) / (-log(freq)+log(n));
			if (v >= threshold && v != 0.) {
				npmi.set(col, v);
			}
		}
		npmi.cr();
	}
	return npmi;
}
