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
#include <phat/stack_access/representations/abstract_pivot.h>

namespace phat { namespace auto_reducing { namespace representations {
    template< typename PivotColumn >
    class abstract_pivot : public stack_access::representations::abstract_pivot< PivotColumn >  {

    protected: 
        std::vector< index > lowest_one_lookup;

    public:
        void _init( index total_num_cols ) {
            stack_access::representations::abstract_pivot< PivotColumn >::_init( total_num_cols );
            lowest_one_lookup.resize( total_num_cols, -1 );
        }

        void _push_col( const column& col, dimension dim ) {
            
            num_cols++;
            const index cur_col = _get_num_cols() - 1;
            column temp_col;
            if( lowest_one_lookup[ cur_col ] == -1 ) {
                pivot_column.clear();
                for( index idx = 0; idx < (index)col.size(); idx++ )
                    pivot_column.add_index( col[ idx ] );
                
                index lowest_one = pivot_column.get_max_index();
                while( lowest_one != -1 && lowest_one_lookup[ lowest_one ] != -1 ) {
                    index column_to_add = lowest_one_lookup[ lowest_one ];
                    for( index entry_idx = offsets[ column_to_add ]; entry_idx < (index)offsets[ column_to_add + 1 ]; entry_idx++ )
                        pivot_column.add_index( entries[ entry_idx ] );
                    lowest_one = pivot_column.get_max_index();
                }
                if( lowest_one != -1 )
                    lowest_one_lookup[ lowest_one ] = cur_col;

                // dump pivot_col into temp_col
                pivot_column.get_col_and_clear( temp_col );
            }

            // dump temp_col into matrix
            dims[ cur_col ] = dim;
            offsets[ cur_col ] = entries.size();
            for( index idx = 0; idx < (index)temp_col.size(); idx++ )
                entries.push_back( temp_col[ idx ] );
            offsets[ cur_col + 1 ] = entries.size();
        }
    };
} } }
