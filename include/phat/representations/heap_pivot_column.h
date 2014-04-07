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
#include <phat/representations/abstract_pivot_column.h>

namespace phat {
    class heap_column {

    protected:
        std::priority_queue< index > data;

        void add_index( const index idx ) {
            data.push( idx );
        }

    public:
        void init( const index total_size ) {
            clear();
        }

        void add_col( const column& col ) {
            for( index idx = 0; idx < (index) col.size(); idx++ )
                add_index( col[ idx ] );
        }

        index get_max_index() {
            if( data.empty() )
                return -1;
            else {
                index max_element = data.top();
                data.pop();
                while( !data.empty() && data.top() == max_element ) {
                    data.pop();
                    if( data.empty() )
                        return -1;
                    else {
                        max_element = data.top();
                        data.pop();
                    }
                }

                data.push( max_element );
                return max_element;
            }
        }

        void get_col_and_clear( column& col ) {
            col.clear();
            while( !data.empty() ) {
                index max_element = data.top();
                data.pop();
                if( data.empty() )
                    col.push_back( max_element );
                else {
                    if( data.top() != max_element )
                        col.push_back( max_element );
                    else
                        data.pop();
                }
            }
            std::reverse( col.begin(), col.end() );
        }

        bool is_empty() {
            return get_max_index() == -1;
        }

		void clear() {
			while( !data.empty() )
                data.pop();
		}

		void remove_max() {
            add_index( get_max_index() );
        }

        void set_col( const column& col  ) {
            clear();
            add_col( col );
        }

        void get_col( column& col  ) {
            get_col_and_clear( col );
            add_col( col );
        }
    };

    typedef abstract_pivot_column< heap_column > heap_pivot_column;
}
