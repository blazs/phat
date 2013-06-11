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
#include <phat/common/read_only_boundary_matrix.h>
#include <phat/random_access/representations/bit_tree_pivot.h>

// interface class for the main data structure -- implementations of the interface can be found in ./representations
namespace phat {
    template< class Representation = bit_tree_pivot_column >
    class random_access_boundary_matrix : public common::read_only_boundary_matrix< Representation >
    {

    // interface functions -- actual implementation and complexity depends on chosen @Representation template
    public:
        // set column @idx to the values contained in @col
        void set_col( index idx, const column& col  ) { rep._set_col( idx, col ); }

        // removes maximal index from given column 
        void remove_max( index idx ) { rep._remove_max( idx ); }

        // adds column @source to column @target'
        void add_to( index source, index target ) { rep._add_to( source, target ); }

        // set dimension of given index
        void set_dim( index idx, dimension dim ) { rep._set_dim( idx, dim ); }

        // clears given column
        void clear( index idx ) { rep._clear( idx ); }

    // operators / constructors
    public:
        random_access_boundary_matrix() {};

        template< class OtherBoundaryMatrix >
        random_access_boundary_matrix( const OtherBoundaryMatrix& other ) {
            *this = other;
        }

        template< typename OtherRepresentation >
        random_access_boundary_matrix< Representation >& operator=( const read_only_boundary_matrix< OtherRepresentation >& other )
        {
            const index nr_of_columns = other.get_num_cols();
            init( nr_of_columns );
            column temp_col;
            for( index cur_col = 0; cur_col <  nr_of_columns; cur_col++ ) {
                this->set_dim( cur_col, other.get_dim( cur_col ) );
                other.get_col( cur_col, temp_col );
                this->set_col( cur_col, temp_col );
            }

            // by convention, always return *this
            return *this;
        }

    // I/O -- independent of chosen 'Representation'
    public:

        // initializes boundary_matrix from (vector<vector>, vector) pair -- untested
        template< typename index_type, typename dimemsion_type >
        void load_vector_vector( const std::vector< std::vector< index_type > >& input_matrix, const std::vector< dimemsion_type >& input_dims ) { 
            const index nr_of_columns = (index)input_matrix.size();
            init( nr_of_columns );
            column temp_col;
            #pragma omp parallel for private( temp_col )
            for( index cur_col = 0; cur_col <  nr_of_columns; cur_col++ ) {
                this->set_dim( cur_col, (dimension)input_dims[ cur_col ] );
                
                index num_rows = input_matrix[ cur_col ].size();
                temp_col.resize( num_rows );
                for( index cur_row = 0; cur_row <  num_rows; cur_row++ )
                    temp_col[ cur_row ] = (index)input_matrix[ cur_col ][ cur_row ];
                this->set_col( cur_col, temp_col );
            }
        }

        // Loads the boundary_matrix from given file in ascii format 
        // Format: each line represents a column, first number is dimension, other numbers are the content of the column.
        // Ignores empty lines and lines starting with a '#'.
        bool load_ascii( std::string filename ) { 
            // first count number of columns:
            std::string cur_line;
            std::ifstream dummy( filename .c_str() );
            if( dummy.fail() )
                return false;

            index number_of_columns = 0;
            while( getline( dummy, cur_line ) ) {
                cur_line.erase(cur_line.find_last_not_of(" \t\n\r\f\v") + 1);
                if( cur_line != "" && cur_line[ 0 ] != '#' )
                    number_of_columns++;

            }
            init( number_of_columns );
            dummy.close();

            std::ifstream input_stream( filename.c_str() );
            if( input_stream.fail() )
                return false;
            
            column temp_col;
            index cur_col = -1;
            while( getline( input_stream, cur_line ) ) {
                cur_line.erase(cur_line.find_last_not_of(" \t\n\r\f\v") + 1);
                if( cur_line != "" && cur_line[ 0 ] != '#' ) {
                    cur_col++;
                    std::stringstream ss( cur_line );
                    
                    int64_t temp_dim;
                    ss >> temp_dim;
                    this->set_dim( cur_col, (dimension) temp_dim );

                    int64_t temp_index;
                    temp_col.clear();
                    while( ss.good() ) {
                        ss >> temp_index;
                        temp_col.push_back( (index)temp_index );
                    }
                    std::sort( temp_col.begin(), temp_col.end() );
                    this->set_col( cur_col, temp_col );
                }
            }

            input_stream.close();
            return true;
        }

        // Loads boundary_matrix from given file 
        // Format: nr_columns % dim1 % N1 % row1 row2 % ...% rowN1 % dim2 % N2 % ...
        bool load_binary( std::string filename ) { 
            std::ifstream input_stream( filename.c_str(), std::ios_base::binary | std::ios_base::in );
            if( input_stream.fail() )
                return false;

            int64_t nr_columns;
            input_stream.read( (char*)&nr_columns, sizeof( int64_t ) );
            init( (index)nr_columns );

            column temp_col;
            for( index cur_col = 0; cur_col < nr_columns; cur_col++ ) {
                int64_t cur_dim;
                input_stream.read( (char*)&cur_dim, sizeof( int64_t ) );
                this->set_dim( cur_col, (dimension) cur_dim );
                int64_t nr_rows;
                input_stream.read( (char*)&nr_rows, sizeof( int64_t ) );
                temp_col.resize( (std::size_t)nr_rows );
                for( index idx = 0; idx < nr_rows; idx++ ) {
                    int64_t cur_row;
                    input_stream.read( (char*)&cur_row, sizeof( int64_t ) );
                    temp_col[ idx ] = (index)cur_row;
                }
                this->set_col( cur_col, temp_col );
            }

            input_stream.close();
            return true;
        }


    };
}
