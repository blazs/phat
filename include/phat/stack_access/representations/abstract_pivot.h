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

#include <phat/common/basic_types.h>
#include <phat/random_access/representations/bit_tree_pivot.h>

namespace phat { namespace stack_access { namespace representations {
    template< typename PivotColumn >
    class abstract_pivot {

    protected:
        std::vector< dimension > dims;
        std::vector< index > entries;
        std::vector< index > offsets;

        mutable bool is_top_column_pivot;
        mutable PivotColumn pivot_column;
        index num_cols;

        void release_pivot_col() {
            column temp_col;
            pivot_column.get_col_and_clear( temp_col );
            entries.resize( offsets[ _get_num_cols() - 1 ] );
            for( index idx = 0; idx < (index)temp_col.size(); idx++ )
                entries.push_back( temp_col[ idx ] );
            offsets[ _get_num_cols() ] = entries.size();
            is_top_column_pivot = false;
        }
        
        void make_pivot_col() {
            pivot_column.clear();
            for( index entry_idx = offsets[ _get_num_cols() - 1 ]; entry_idx < (index)offsets[ _get_num_cols() ]; entry_idx++ )
                pivot_column.add_index( entries[ entry_idx ] );
            is_top_column_pivot = true;
        }


    public:
        // overall number of cells in boundary_matrix
        index _get_num_cols() const {
            return num_cols; 
        }

        void _init( index total_num_cols ) {
            offsets.resize( total_num_cols + 1, -1 );
            dims.resize( total_num_cols, -1 );
            pivot_column.init( total_num_cols );
            is_top_column_pivot = false;
            num_cols = 0;
        }

        // dimension of given index
        dimension _get_dim( index idx ) const { 
            return dims[ idx ]; 
        }

        // replaces(!) content of 'col' with boundary of given index
        void _get_col( index idx, column& col  ) const {
            if( idx == _get_num_cols() - 1 && is_top_column_pivot ) {
                pivot_column.get_col( col );
            } else {
                col.clear();
                col.reserve( offsets[ idx + 1 ] - offsets[ idx ] );
                for( index entry_idx = offsets[ idx ]; entry_idx < (index)offsets[ idx + 1 ]; entry_idx++ )
                    col.push_back( entries[ entry_idx ] );
            }
        }

        // true iff boundary of given idx is empty
        bool _is_empty( index idx ) const { 
            if( idx == _get_num_cols() - 1 && is_top_column_pivot )
                return pivot_column.is_empty();
            else 
                return offsets[ idx ] == offsets[ idx + 1 ];
        }

        // largest row index of given column idx (new name for lowestOne())
        index _get_max_index( index idx ) const {
            if( idx == _get_num_cols() - 1 && is_top_column_pivot )
                return pivot_column.get_max_index();
            else 
                return (offsets[ idx ] == offsets[ idx + 1 ]) ? -1 : entries[ offsets[ idx + 1 ] - 1 ];
        }

        // syncronizes all data structures (essential for openmp stuff)
        void _sync() {}

        void _push_col( const column& col, dimension dim  ) {
            if( is_top_column_pivot )
                release_pivot_col();

            dims[ _get_num_cols() ] = dim;
            
            offsets[ _get_num_cols() ] = entries.size();
            for( index idx = 0; idx < (index)col.size(); idx++ )
                entries.push_back( col[ idx ] );
            num_cols++;
            offsets[ _get_num_cols() ] =  entries.size();
            
        }

        void _add_to_top( index source ) {
            if( !is_top_column_pivot )
                make_pivot_col();

            for( index entry_idx = offsets[ source ]; entry_idx < (index)offsets[ source + 1 ]; entry_idx++ )
                pivot_column.add_index( entries[ entry_idx ] );
        }
    };
} } }
