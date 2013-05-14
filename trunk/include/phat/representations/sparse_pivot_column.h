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
    class sparse_column {

    protected:
        std::set< index > m_data;

        void add_index( const index idx ) {
            std::pair< std::set< index >::iterator, bool > result = m_data.insert( idx );
            if( result.second == false )
                m_data.erase( result.first );
        }

    public:
        void init( const index total_size ) {
            m_data.clear(); 
        }

        void add_col( const column& col ) {
            for( index idx = 0; idx < (index) col.size(); idx++ )
                add_index( col[ idx ] );
        }

        index get_max_index() {
            return m_data.empty() ? -1 : *m_data.rbegin();
        }

        void get_col_and_clear( column& col ) {
            col.clear();
            col.assign( m_data.begin(), m_data.end() );
            m_data.clear();
        }

        bool is_empty() {
            return m_data.empty();
        }

		void clear() {
			m_data.clear();
		}

		void remove_max() {
            add_index( get_max_index() );
        }

        void set_col( const column& col  ) {
            clear();
            add_col( col );
        }

        void get_col( column& col  ) {
            col.clear();
            get_col_and_clear( col );
            add_col( col );
        }
    };

    typedef abstract_pivot_column< sparse_column > sparse_pivot_column;
}
