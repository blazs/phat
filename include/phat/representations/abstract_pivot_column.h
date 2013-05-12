/*  Copyright 2013 IST Austria
    Contributed by: Ulrich Bauer, Michael Kerber, Jan Reininghaus

    This file is part of PHAT.

    PHAT is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PHAT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with PHAT.  If not, see <http://www.gnu.org/licenses/>. */

#pragma once

#include <phat/helpers/misc.h>
#include <phat/representations/vector_vector.h>

namespace phat {

    // Note: We could even make the rep generic in the underlying Const representation
    //       But I cannot imagine that anything else than vector<vector<index>> would
    //       make sense
    template< typename PivotColumn >
    class abstract_pivot_column : public vector_vector {
    public:
        
    protected:
        typedef vector_vector Base;
        typedef PivotColumn pivot_column;

        // For parallization purposes, it could be more than one full column
        mutable thread_local_storage< pivot_column > _pivot_columns;
        mutable thread_local_storage< index > pos_of_pivot_columns;

        pivot_column& get_pivot_column() const {
            return _pivot_columns();
        }

        bool is_represented_by_pivot_column( index idx ) const {
            return pos_of_pivot_columns() == idx;
        }
      
        void unset_pos_of_pivot_column() {
            index idx = pos_of_pivot_columns();
            if( idx != -1 ) {
                _pivot_columns().get_column_and_clear( this->matrix[ idx ] );
            }
            pos_of_pivot_columns() = -1;
        }
        
        void represent_by_pivot_column( index idx ) {
            pos_of_pivot_columns() = idx;
            get_pivot_column().add_column( matrix[ idx ] );
        }

    public:  

        void _set_num_cols( index nr_of_columns ) {
            #pragma omp parallel for
            for( int tid = 0; tid < omp_get_num_threads(); tid++ ) {
                _pivot_columns[ tid ].init( nr_of_columns );
                pos_of_pivot_columns[ tid ] = -1;
            }
            Base::_set_num_cols( nr_of_columns );
        }
        // replaces(!) content of 'col' with boundary of given index
        void _get_col( index idx, column& col  ) const {
            is_represented_by_pivot_column( idx ) ? get_pivot_column().get_col( col ) : Base::_get_col( idx, col );
        }
        
        // true iff boundary of given idx is empty
        bool _is_empty( index idx ) const {
            return is_represented_by_pivot_column( idx ) ? get_pivot_column().empty() : Base::_is_empty( idx );
        }

        // largest row index of given column idx (new name for lowestOne())
        index _get_max_index( index idx ) const {
			return is_represented_by_pivot_column( idx ) ? get_pivot_column().max_index() : Base::_get_max_index( idx );
        }

        // adds column 'source' to column 'target'
        void _add_to( index source, index target ) {
            if( !is_represented_by_pivot_column( target ) ) {
                unset_pos_of_pivot_column();
                represent_by_pivot_column( target );
            } 
            get_pivot_column().add_column( matrix[source] );
        }

        // clears given column
        void _clear( index idx ) {
			is_represented_by_pivot_column( idx ) ? get_pivot_column().clear() : Base::_clear( idx );
        }

        void _set_col( index idx, const column& col  ) {
            is_represented_by_pivot_column( idx ) ? get_pivot_column().set_col( col ) : Base::_set_col( idx, col );
        }

        // removes the maximal index of a column
        void _remove_max( index idx ) {
			is_represented_by_pivot_column( idx ) ? get_pivot_column().remove_max() : Base::_remove_max( idx );
        }

        // syncronizes all data structures (essential for openmp stuff)
        // has to be called before and after any multithreaded access!
        void _sync() { 
            #pragma omp parallel for
            for( int tid = 0; tid < omp_get_num_threads(); tid++ )
                unset_pos_of_pivot_column();
        } 

    };
}


