#ifndef SPMTRX_H
#define SPMTRX_H

namespace spmtrx {
	class spmtrx {
		public:
			virtual spmtrx& operator=(spmtrx&& a) noexcept;
			virtual spmtrx& operator=(spmtrx& a) = 0;
			virtual spmtrx operator+(const spmtrx& a) = 0;
			virtual spmtrx operator+(const spmtrx&& a) = 0;
			virtual spmtrx operator-(const spmtrx& a) = 0;
			virtual spmtrx operator-(const spmtrx&& a) = 0;
			virtual spmtrx& operator+=(const spmtrx& a) = 0;
			virtual spmtrx& operator+=(const spmtrx&& a) = 0;
			virtual spmtrx& operator-=(const spmtrx& a) = 0;
			virtual spmtrx& operator-=(const spmtrx&& a) = 0;
			virtual spmtrx operator*(const spmtrx& a) = 0;
			virtual spmtrx operator*(const spmtrx&& a) = 0;
			virtual spmtrx operator*(double c) = 0;
			virtual spmtrx& operator*=(const spmtrx& a) = 0;
			virtual spmtrx& operator*=(const spmtrx&& a) = 0;
			virtual spmtrx& operator*=(double c) = 0;
			virtual double get(int row, int col) = 0;
			virtual int get_head(int i) = 0;
			virtual double get_val(int id) = 0;
			virtual int get_row_size() = 0;
			virtual int get_col_size() = 0;
			virtual spmtrx transpose() = 0;
			virtual spmtrx diag() = 0;
			virtual spmtrx inv_diag() = 0;
			virtual spmtrx sqrt_inv_diag() = 0;
			virtual spmtrx diag_adjacency() = 0;
			virtual spmtrx inv_diag_adjacency() = 0;
			virtual spmtrx sqrt_inv_diag_adjacency() = 0;
			virtual spmtrx npmi(double threshold) = 0;
			virtual double trace() = 0;
			virtual double self_dot_prod(int i, int j) = 0;
			virtual void add(int head_id, double v) = 0;
			virtual void cr() = 0;
			virtual void dump() = 0;
			virtual const_iterator begin() = 0;
			virtual const_iterator end() = 0;
		protected:
			friend class boost::serialization::access;
			template<class Archive>
			virtual void serialize(Archive& ar, unsigned int version) = 0;
			
	};
}
#endif /** SPMTRX_H **/
