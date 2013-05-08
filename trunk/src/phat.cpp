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

#include <phat/compute_persistence_pairs.h>

#include <phat/representations/vector_vector.h>
#include <phat/representations/vector_set.h>
#include <phat/representations/vector_list.h>
#include <phat/representations/sparse_pivot_column.h>
#include <phat/representations/full_pivot_column.h>
#include <phat/representations/bit_tree_pivot_column.h>

#include <phat/algorithms/twist_reduction.h>
#include <phat/algorithms/standard_reduction.h>
#include <phat/algorithms/row_reduction.h>
#include <phat/algorithms/chunk_reduction.h>


enum Representation_type  {VECTOR_VECTOR, VECTOR_SET, SPARSE_PIVOT_COLUMN, FULL_PIVOT_COLUMN, BIT_TREE_PIVOT_COLUMN, VECTOR_LIST};
enum Algorithm_type  {STANDARD, TWIST, ROW, CHUNK, CHUNK_SEQUENTIAL };

void print_help() {
    std::cerr << "Usage: " << "phat " << "[options] input_filename output_filename" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << std::endl;
    std::cerr << "--ascii   --  use ascii file format" << std::endl;
    std::cerr << "--binary  --  use binary file format (default)" << std::endl;
    std::cerr << "--help    --  prints this screen" << std::endl;
    std::cerr << "--verbose --  verbose output" << std::endl;
    std::cerr << "--dualize   --  use dualization approach" << std::endl;
    std::cerr << "--vector_vector, --vector_set, --vector_list, --full_pivot_column, --sparse_pivot_column, --bit_tree_pivot_column  --  selects a representation data structure for boundary matrices (default is '--bit_tree_pivot_column')" << std::endl;
    std::cerr << "--standard, --twist, --chunk, --chunk_sequential, --row  --  selects a reduction algorithm (default is '--twist')" << std::endl;
}

void print_help_and_exit() {
    print_help();
    exit( EXIT_FAILURE );
}

void parse_command_line( int argc, char** argv, bool& use_binary, Representation_type& representation, Algorithm_type& algorithm,
                         std::string& input_filename, std::string& output_filename, bool& verbose, bool& dualize) {

    if( argc < 3 ) print_help_and_exit();

    input_filename = argv[ argc - 2 ];
    output_filename = argv[ argc - 1 ];

    for( int idx = 1; idx < argc - 2; idx++ ) {
        const std::string option = argv[ idx ];

        if( option == "--ascii" ) use_binary = false;
        else if( option == "--binary" ) use_binary = true;
        else if( option == "--dualize" ) dualize = true;
        else if( option == "--vector_vector" ) representation = VECTOR_VECTOR;
        else if( option == "--vector_set" ) representation = VECTOR_SET;
        else if( option == "--vector_list" ) representation = VECTOR_LIST;
        else if( option == "--full_pivot_column" )  representation = FULL_PIVOT_COLUMN;
        else if( option == "--bit_tree_pivot_column" )  representation = BIT_TREE_PIVOT_COLUMN;
        else if( option == "--sparse_pivot_column" ) representation = SPARSE_PIVOT_COLUMN;
        else if( option == "--standard" ) algorithm = STANDARD;
        else if( option == "--twist" ) algorithm = TWIST;
        else if( option == "--row" ) algorithm = ROW;
        else if( option == "--chunk" ) algorithm = CHUNK;
        else if( option == "--chunk_sequential" ) algorithm = CHUNK_SEQUENTIAL;
        else if( option == "--verbose" ) verbose = true;
        else if( option == "--help" ) print_help_and_exit();
        else print_help_and_exit();
    }
}

#define LOG(msg) if( verbose ) std::cout << msg << std::endl;

template<typename Representation, typename Algorithm>
void compute_pairing( std::string input_filename, std::string output_filename, bool use_binary, bool verbose, bool dualize ) {

    phat::boundary_matrix< Representation > matrix;
    bool read_successful;

    double read_timer = omp_get_wtime();
    if( use_binary ) {
        LOG( "Reading input file " << input_filename << " in binary mode" )
        read_successful = matrix.load_binary( input_filename );
    } else {
        LOG( "Reading input file " << input_filename << " in ascii mode" )
        read_successful = matrix.load_ascii( input_filename );
    }
    LOG( "Reading input file took " << omp_get_wtime() - read_timer <<"s" )

    if( !read_successful ) {
        std::cerr << "Error opening file " << input_filename << std::endl;
        print_help_and_exit();
    }

    double pairs_timer = omp_get_wtime();
    phat::persistence_pairs pairs;
    LOG( "Computing persistence pairs ..." )
    if( dualize )
        phat::compute_persistence_pairs_dualized< Algorithm > ( pairs, matrix );
    else
        phat::compute_persistence_pairs < Algorithm > ( pairs, matrix );
    LOG( "Computing persistence pairs took " << omp_get_wtime() - pairs_timer <<"s" )

    double write_timer = omp_get_wtime();
    if( use_binary ) {
        LOG( "Writing output file " << output_filename << " in binary mode ..." )
        pairs.save_binary( output_filename );
    } else {
        LOG( "Writing output file " << output_filename << " in ascii mode ..." )
        pairs.save_ascii( output_filename );
    }
    LOG( "Writing output file took " << omp_get_wtime() - write_timer <<"s" )
}

#define COMPUTE_PAIRING(Representation) \
    switch( algorithm ) { \
    case STANDARD: compute_pairing< phat::Representation, phat::standard_reduction> ( input_filename, output_filename, use_binary, verbose, dualize ); break; \
    case TWIST: compute_pairing< phat::Representation, phat::twist_reduction> ( input_filename, output_filename, use_binary, verbose, dualize ); break; \
    case ROW: compute_pairing< phat::Representation, phat::row_reduction >( input_filename, output_filename, use_binary, verbose, dualize ); break; \
    case CHUNK: compute_pairing< phat::Representation, phat::chunk_reduction >( input_filename, output_filename, use_binary, verbose, dualize ); break; \
    case CHUNK_SEQUENTIAL: int num_threads = omp_get_max_threads(); \
                           omp_set_num_threads( 1 ); \
                           compute_pairing< phat::Representation, phat::chunk_reduction >( input_filename, output_filename, use_binary, verbose, dualize ); break; \
                           omp_set_num_threads( num_threads ); \
                           break; \
    }

int main( int argc, char** argv )
{
    bool use_binary = true; // interpret input as binary or ascii file
    Representation_type representation = BIT_TREE_PIVOT_COLUMN; // representation class
    Algorithm_type algorithm = TWIST; // reduction algorithm
    std::string input_filename; // name of file that contains the boundary matrix
    std::string output_filename; // name of file that will contain the persistence pairs
    bool verbose = false; // print timings / info
    bool dualize = false; // toggle for dualization approach

    parse_command_line( argc, argv, use_binary, representation, algorithm, input_filename, output_filename, verbose, dualize );

    switch( representation ) {
    case VECTOR_VECTOR: COMPUTE_PAIRING(vector_vector) break;
    case VECTOR_SET: COMPUTE_PAIRING(vector_set) break;
    case VECTOR_LIST: COMPUTE_PAIRING(vector_list) break;
    case FULL_PIVOT_COLUMN: COMPUTE_PAIRING(full_pivot_column) break;
    case BIT_TREE_PIVOT_COLUMN: COMPUTE_PAIRING(bit_tree_pivot_column) break;
    case SPARSE_PIVOT_COLUMN: COMPUTE_PAIRING(sparse_pivot_column) break;
    }
}
