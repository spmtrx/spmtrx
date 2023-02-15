#include"spmtrx/spmcrs.h"
#include<cstdlib>

using namespace std;
using namespace spm;

spmcrs::spmcrs(): _data(new vector<double>), _row(new vector<int>(1,0)), _col(new vector<int>), _colsize(0) {
}

spmcrs::spmcrs(spmcrs& a): _data(new vector<double>), _row(new vector<int>(1,0)), _col(new vector<int>), _colsize(0) {
	for (auto it = a.begin(); it != a.end(); ++it) {
		for (auto j = *it; j < *(it+1); ++j) {
			set((*a._col)[j], (*a._data)[j]);
		}
		cr();
	}
}

spmcrs::spmcrs(spmcrs&& a) noexcept {
	_data = a._data;
	_row = a._row;
	_col = a._col;
	_colsize = a._colsize;
	a._data = nullptr;
	a._row = nullptr;
	a._col = nullptr;
}

spmcrs::~spmcrs() {
}

void spmcrs::set(unsigned int col, double v) {
	_col->push_back(col);
	_data->push_back(v);
	_colsize = max(_colsize, (int)col+1);
}

void spmcrs::cr() {
	_row->push_back(_data->size());
}

void spmcrs::dump() {
	auto i = 0;
	for (auto it = begin(); it != end(); ++it, ++i) {
		for (auto j = *it; j != *(it+1); ++j) {
			cout << i << "," << (*_col)[j] << ":" << (*_data)[j] << "\t";
		}
		if (*it != *(it+1))
			cout << endl;
	}
}

vector<int>::iterator spmcrs::begin() {
	return _row->begin();
}

vector<int>::iterator spmcrs::end() {
	return _row->end()-1;
}

vector<int>::const_iterator spmcrs::cbegin() {
	return _row->cbegin();
}

vector<int>::const_iterator spmcrs::cend() {
	return _row->cend()-1;
}

unsigned int spmcrs::get_col(unsigned int id) {
	if (id >= _col->size())
		return 0;
	return (*_col)[id];
}

double spmcrs::get_val(unsigned int id) {
	if (id >= _data->size())
		return 0;
	return (*_data)[id];
}

double spmcrs::get_val(unsigned int row, unsigned int col) {
	if (row >= _row->size()-1)
		return 0;
	auto head = *(_row->begin()+row);
	auto tail = *(_row->begin()+row+1);
	int c = 0;
	while (head != tail) {
		c = (head + tail)/2;
		if (head == c && (*_col)[c] != col) {
			return 0; // not found
		} else if ((*_col)[c] == col) {
			return (*_data)[c];
		} else if ((*_col)[c] > col) {
			tail = c;
		} else if ((*_col)[c] < col) {
			head = c;
		}
	}
	return 0;
}

int spmcrs::get_row_size() {
	return _row->size()-1;
}

int spmcrs::get_col_size() {
	return _colsize;
}

spmcrs spmcrs::ones(unsigned int size) {
	spmcrs unit;
	unit._data->resize(size, 1);
	unit._row->resize(size+1, 0);
	unit._col->resize(size, 0);
	unit._colsize = 1;
	for (auto i = 0; i < size+1; ++i)
		(*unit._row)[i] = i;
	return unit;
}

spmcrs spmcrs::transpose() {
	spmcrs t;
	t._data->resize(_data->size());
	t._col->resize(_data->size());
	t._row->resize(_colsize+1);
	t._colsize = _row->size();

	vector<int> t_non_zero_in_row(_colsize+1,0);
	for (auto i = 0; i < _data->size(); ++i) {
		t_non_zero_in_row[(*_col)[i]]++;
	}
	for (auto i = 0, n = 0; i < _colsize+1; ++i) {
		(*t._row)[i] = n;
		n += t_non_zero_in_row[i];
	}
	vector<int> count(_colsize,0);
	auto it = begin();
	for (auto i = 0; it != end(); ++it, ++i) {
		for (auto j = *it; j < *(it+1); ++j) {
			int t_id = (*t._row)[(*_col)[j]]+count[(*_col)[j]]++;
			(*t._data)[t_id] = (*_data)[j];
			(*t._col)[t_id] = i;
		}
	}
	return t;
}

spmcrs spmcrs::diag() {
	spmcrs d;
	auto it = begin();
	for (auto i = 0; it != end(); ++i, ++it) {
		double v = 0;
		for (auto j = *it; j < *(it+1); ++j) {
			v += (*_data)[j];
		}
		if (v != 0.)
			d.set(i, v);
		d.cr();
	}
	return d;
}

spmcrs spmcrs::inv_diag() {
	spmcrs id = diag();
	for (auto it = id._data->begin(); it != id._data->end(); ++it) {
		*it = 1. / *it;
	}
	return id;
}

spmcrs spmcrs::sqrt_inv_diag() {
	spmcrs sid = diag();
	for (auto it = sid._data->begin(); it != sid._data->end(); ++it) {
		*it = 1./sqrt(*it);
	}
	return sid;
}

spmcrs spmcrs::diag_adjacency() {
	spmcrs t = transpose();
	spmcrs e = ones(t._colsize);
	spmcrs a = *this*(t*e);
	return a.diag();
}

spmcrs spmcrs::inv_diag_adjacency() {
	spmcrs t = transpose();
	spmcrs e = ones(t._colsize);
	spmcrs a = *this*(t*e);
	return a.inv_diag();
}

spmcrs spmcrs::sqrt_inv_diag_adjacency() {
	spmcrs t = transpose();
	spmcrs e = ones(t._colsize);
	spmcrs a = *this*(t*e);
	return a.sqrt_inv_diag();
}

double spmcrs::trace() {
	double tr = 0;
	auto i = 0;
	for  (auto it = begin(); it != end(); ++it, ++i)
		tr += get_val(i,i);
	return tr;
}

double spmcrs::self_dot_prod(unsigned int i, unsigned int j) {
	if (max(i, j) >= get_row_size())
		return 0;
	double v = 0;
	auto it_i = begin()+i;
	auto it_j = begin()+j;
	auto i_id = *it_i;
	auto j_id = *it_j;
	while(i_id < *(it_i+1) && j_id < *(it_j+1)) {
		auto icol = (*_col)[i_id];
		auto jcol = (*_col)[j_id];
		if (icol < jcol) {
			++i_id;
		} else if (icol > jcol) {
			++j_id;
		} else if (icol == jcol) {
			v += (*_data)[i_id]*(*_data)[j_id];
		}
	}
	return v;
}



spmcrs& spmcrs::operator=(spmcrs& x) {
	_data = x._data;
	_row = x._row;
	_col = x._col;
	_colsize = x._colsize;
	return *this;
}

spmcrs& spmcrs::operator=(spmcrs&& x) noexcept {
	_data = x._data;
	_row = x._row;
	_col = x._col;
	_colsize = x._colsize;
	x._data = nullptr;
	x._row = nullptr;
	x._col = nullptr;
	return *this;
}

spmcrs spmcrs::operator+(spmcrs& b) {
	spmcrs y;
	auto a_it = begin();
	auto b_it = b.begin();
	while (a_it != end() && b_it != b.end()) {
		for (auto a_j = *a_it, b_j = *b_it; (a_j < *(a_it+1)) || (b_j < *(b_it+1));) {
			if (a_j == *(a_it+1) && b_j < *(b_it+1)) {
				y.set((*b._col)[b_j], (*b._data)[b_j]);
				++b_j;
			} else if (a_j < *(a_it+1) && b_j == *(b_it+1)) {
				y.set((*_col)[a_j], (*_data)[a_j]);
				++a_j;
			} else if ((*_col)[a_j] > (*b._col)[b_j]) {
				y.set((*b._col)[b_j], (*b._data)[b_j]);
				++b_j;
			} else if ((*_col)[a_j] < (*b._col)[b_j]) {
				y.set((*_col)[a_j], (*_data)[a_j]);
				++a_j;

			} else if ((*_col)[a_j] == (*b._col)[b_j]) {
				double v = (*_data)[a_j] + (*b._data)[b_j];
				if (v != 0.)
					y.set((*_col)[a_j], v);
				++a_j, ++b_j;
			}
		}
		y.cr();
		++a_it, ++b_it;
	}
	for (; a_it != cend(); ++a_it) {
		for (auto j = *a_it; j < *(a_it+1); ++j)
			y.set((*_col)[j], (*_data)[j]);
		y.cr();
	}
	for (; b_it != b.cend(); ++b_it) {
		for (auto j = *b_it; j < *(b_it+1); ++j)
			y.set((*b._col)[j], (*b._data)[j]);
		y.cr();
	}
	return y;
}

spmcrs spmcrs::operator+(spmcrs&& b) {
	spmcrs y = *this+b;
	return y;
}

spmcrs spmcrs::operator-(spmcrs& b) {
	spmcrs y;
	auto a_it = begin();
	auto b_it = b.begin();
	while (a_it != end() && b_it != b.end()) {
		for (auto a_j = *a_it, b_j = *b_it; (a_j < *(a_it+1)) || (b_j < *(b_it+1));) {
			if (a_j == *(a_it+1) && b_j < *(b_it+1)) {
				y.set((*b._col)[b_j], -(*b._data)[b_j]);
				++b_j;
			} else if (a_j < *(a_it+1) && b_j == *(b_it+1)) {
				y.set((*_col)[a_j], (*_data)[a_j]);
				++a_j;
			} else if ((*_col)[a_j] > (*b._col)[b_j]) {
				y.set((*b._col)[b_j], -(*b._data)[b_j]);
				++b_j;
			} else if ((*_col)[a_j] < (*b._col)[b_j]) {
				y.set((*_col)[a_j], (*_data)[a_j]);
				++a_j;

			} else if ((*_col)[a_j] == (*b._col)[b_j]) {
				double v = (*_data)[a_j] - (*b._data)[b_j];
				if (v != 0.)
					y.set((*_col)[a_j], v);
				++a_j, ++b_j;
			}
		}
		y.cr();
		++a_it, ++b_it;
	}
	for (; a_it != end(); ++a_it) {
		for (auto j = *a_it; j < *(a_it+1); ++j)
			y.set((*_col)[j], (*_data)[j]);
		y.cr();
	}
	for (; b_it != b.end(); ++b_it) {
		for (auto j = *b_it; j < *(b_it+1); ++j)
			y.set((*b._col)[j], (*b._data)[j]);
		y.cr();
	}
	return y;
}

spmcrs spmcrs::operator-(spmcrs&& b) {
	spmcrs y = *this-b;
	return y;
}

spmcrs& spmcrs::operator+=(spmcrs& a) {
	*this = *this+a;
	return *this;
}

spmcrs& spmcrs::operator+=(spmcrs&& a) {
	*this = *this+a;
	return *this;
}

spmcrs& spmcrs::operator-=(spmcrs& a) {
	*this = *this-a;
	return *this;
}

spmcrs& spmcrs::operator-=(spmcrs&& a) {
	*this = *this-a;
	return *this;
}

spmcrs spmcrs::operator*(double c) {
	spmcrs x(*this);
	for (auto it = x._data->begin(); it != x._data->end(); ++it)
		*it *= c;
	return x;
}

spmcrs& spmcrs::operator*=(double c) {
	for (auto it = _data->begin(); it != _data->end(); ++it)
		*it *= c;
	return *this;
}

spmcrs spmcrs::operator/(double c) {
	return *this*1./c;
}

spmcrs& spmcrs::operator/=(double c) {
	return *this*=1./c;
}

spmcrs spmcrs::operator*(spmcrs& a) {
	spmcrs x;
	auto it = begin();
	for (auto i = 0; it != end(); ++it, ++i) {
		for (auto j = 0; j < a._colsize; ++j) {
			double v = 0;
			for (auto k = *it; k < *(it+1); ++k) {
				v += (*_data)[k] * a.get_val((*_col)[k], j);
			}
			if (v != 0.)
				x.set(j, v);
		}
		x.cr();
	}
	return x;
}

spmcrs spmcrs::operator*(spmcrs&& a) {
	return *this*a;
}

spmcrs& spmcrs::operator*=(spmcrs& a) {
	*this = *this*a;
	return *this;
}

spmcrs& spmcrs::operator*=(spmcrs&& a) {
	*this = *this*a;
	return *this;
}
