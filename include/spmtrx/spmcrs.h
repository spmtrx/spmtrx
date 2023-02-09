#ifndef SPMCRS_H
#define SPMCRS_H

#include<vector>
#include<memory>
#include<boost/serialization/vector.hpp>
#include<boost/archive/text_oarchive.hpp>
#include<boost/archive/text_iarchive.hpp>

namespace spm {
	class spmcrs {
		public:
			spmcrs();
			spmcrs(spmcrs& a);
			spmcrs(spmcrs&& a) noexcept;
			~spmcrs();
			spmcrs& operator=(spmcrs&& a) noexcept;
			spmcrs& operator=(spmcrs& a);
			spmcrs operator+(spmcrs& a);
			spmcrs operator+(spmcrs&& a);
			spmcrs operator-(spmcrs& a);
			spmcrs operator-(spmcrs&& a);
			spmcrs& operator+=(spmcrs& a);
			spmcrs& operator+=(spmcrs&& a);
			spmcrs& operator-=(spmcrs& a);
			spmcrs& operator-=(spmcrs&& a);
			spmcrs operator*(spmcrs& a);
			spmcrs operator*(spmcrs&& a);
			spmcrs operator*(double c);
			spmcrs& operator*=(spmcrs& a);
			spmcrs& operator*=(spmcrs&& a);
			spmcrs& operator*=(double c);
			spmcrs operator/(double c);
			spmcrs& operator/=(double c);
			unsigned int get_col(unsigned int id);
			double get_val(unsigned int id);
			double get_val(unsigned int row, unsigned int col);
			int get_row_size();
			int get_col_size();
			spmcrs ones(unsigned int size);
			spmcrs transpose();
			spmcrs diag();
			spmcrs inv_diag();
			spmcrs sqrt_inv_diag();
			spmcrs diag_adjacency();
			spmcrs inv_diag_adjacency();
			spmcrs sqrt_inv_diag_adjacency();
			spmcrs npmi(double threshold);
			double trace();
			double self_dot_prod(unsigned int i, unsigned int j);
			void set(unsigned int col, double v);
			void cr();
			void dump();
			std::vector<int>::iterator begin();
			std::vector<int>::iterator end();
			std::vector<int>::const_iterator cbegin();
			std::vector<int>::const_iterator cend();
		protected:
			std::shared_ptr<std::vector<double> > _data;
			std::shared_ptr<std::vector<int> > _row;
			std::shared_ptr<std::vector<int> > _col;
			int _colsize;
		private:
			friend class boost::serialization::access;
			template<class Archive>
			void serialize(Archive& ar, unsigned int version) {
				ar & *_data;
				ar & *_row;
				ar & *_col;
				ar & _colsize;
			}
			
	};
}
#endif /** SPMTRX_H **/
