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
#include <phat/stack_access/boundary_matrix.h>

// interface class for the main data structure -- implementations of the interface can be found in ./representations
namespace phat { namespace auto_reducing {
    template< class Representation = representations::bit_tree_pivot >
    class boundary_matrix : public stack_access::boundary_matrix< Representation > {
    public:
        boundary_matrix() {};

        template< class OtherBoundaryMatrix >
        boundary_matrix( const OtherBoundaryMatrix& other ) {
            *this = other;
        }

        template< typename OtherRepresentation >
        boundary_matrix< Representation >& operator=( const common::const_boundary_matrix< OtherRepresentation >& other )
        {
            if( (void*)this != (void*)&other ) {
                const index nr_of_columns = other.get_num_cols();
                init( nr_of_columns );
                column temp_col;
                for( index cur_col = 0; cur_col <  nr_of_columns; cur_col++ ) {
                    other.get_col( cur_col, temp_col );
                    this->push_col( temp_col, other.get_dim( cur_col ) );
                }
            }

            // by convention, always return *this
            return *this;
        }
    };
} }
