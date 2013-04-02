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
#include <phat/representations/sparse_pivot_column.h>
#include <phat/representations/full_pivot_column.h>
#include <phat/representations/bit_tree_pivot_column.h>

#include <phat/algorithms/twist_reduction.h>
#include <phat/algorithms/standard_reduction.h>
#include <phat/algorithms/row_reduction.h>
#include <phat/algorithms/chunk_reduction.h>


enum Representation_type  {VEC_VEC, VEC_SET, SPARSE_PIVOT, FULL_PIVOT, BIT_TREE_PIVOT};
enum Algorithm_type  {STANDARD, TWIST, ROW, CHUNK };

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
    std::cerr << "--vec-vec, --vec-set, --full-pivot, --sparse-pivot, --bit-tree-pivot  --  selects a representation data structure for boundary matrices (default is '--sparse-pivot')" << std::endl;
    std::cerr << "--standard, --twist, --chunk, --row  --  selects a reduction algorithm (default is '--twist')" << std::endl;
}

void print_help_and_exit() {
    print_help();
    exit( EXIT_FAILURE );
}

void parse_command_line( int argc, char** argv, bool& use_binary, Representation_type& rep_type, Algorithm_type& reduction,
                         std::string& input_filename, std::string& output_filename, bool& verbose, bool& dualize ) {

    if( argc < 3 ) print_help_and_exit();

    input_filename = argv[ argc - 2 ];
    output_filename = argv[ argc - 1 ];

    for( int idx = 1; idx < argc - 2; idx++ ) {
        const std::string option = argv[ idx ];

        if( option == "--ascii" ) use_binary = false;
        else if( option == "--binary" ) use_binary = true;
        else if( option == "--dualize" ) dualize = true;
        else if( option == "--vec-vec" ) rep_type = VEC_VEC;
        else if( option == "--vec-set" ) rep_type = VEC_SET;
        else if( option == "--full-pivot" )  rep_type = FULL_PIVOT;
		else if( option == "--bit-tree-pivot" )  rep_type = BIT_TREE_PIVOT;
        else if( option == "--sparse-pivot" ) rep_type = SPARSE_PIVOT;
        else if( option == "--standard" ) reduction = STANDARD;
        else if( option == "--twist" ) reduction = TWIST;
        else if( option == "--row" ) reduction = ROW;
        else if( option == "--chunk" ) reduction = CHUNK;
        else if( option == "--verbose" ) verbose = true;
        else if( option == "--help" ) print_help_and_exit();
        else print_help_and_exit();
    }
}

#define LOG(msg) if( verbose ) std::cout << msg << std::endl;

template<typename Representation, typename Algorithm>
void generic_compute_pairing( std::string input_filename,
                              std::string output_filename,
	                          bool use_binary,
	                          bool verbose,
                              bool dualize ) {

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

#define CALL_GENERIC_CODE(rep,alg) generic_compute_pairing < rep, alg >( input_filename, output_filename, use_binary, verbose, dualize );

int main( int argc, char** argv )
{
    bool use_binary = true; // interpret input as binary or ascii file
    Representation_type rep_type = SPARSE_PIVOT; // representation class
    Algorithm_type reduction = TWIST; // reduction algorithm
    std::string input_filename; // name of file that contains the boundary matrix
    std::string output_filename; // name of file that will contain the persistence pairs
    bool verbose = false; // print timings / info
    bool dualize = false; // toggle for dualization approach

    parse_command_line( argc, argv, use_binary, rep_type, reduction, input_filename, output_filename, verbose, dualize );

    switch( rep_type ) {
    case VEC_VEC:       switch( reduction ) {
                        case STANDARD: CALL_GENERIC_CODE(phat::vector_vector, phat::standard_reduction) break;
                        case TWIST: CALL_GENERIC_CODE(phat::vector_vector, phat::twist_reduction) break;
                        case ROW: CALL_GENERIC_CODE(phat::vector_vector, phat::row_reduction) break;
                        case CHUNK: CALL_GENERIC_CODE(phat::vector_vector, phat::chunk_reduction) break;
                        } break;

    case VEC_SET:       switch( reduction ) {
                        case STANDARD: CALL_GENERIC_CODE(phat::vector_set, phat::standard_reduction) break;
                        case TWIST: CALL_GENERIC_CODE(phat::vector_set, phat::twist_reduction) break;
                        case ROW: CALL_GENERIC_CODE(phat::vector_set, phat::row_reduction) break;
                        case CHUNK: CALL_GENERIC_CODE(phat::vector_set, phat::chunk_reduction) break;
                        } break;

    case FULL_PIVOT:    switch( reduction ) {
                        case STANDARD: CALL_GENERIC_CODE(phat::full_pivot_column, phat::standard_reduction) break;
                        case TWIST: CALL_GENERIC_CODE(phat::full_pivot_column, phat::twist_reduction) break;
                        case ROW: CALL_GENERIC_CODE(phat::full_pivot_column, phat::row_reduction) break;
                        case CHUNK: CALL_GENERIC_CODE(phat::full_pivot_column, phat::chunk_reduction) break;
                        } break;

	case BIT_TREE_PIVOT:  switch( reduction ) {
						case STANDARD: CALL_GENERIC_CODE(phat::bit_tree_pivot_column, phat::standard_reduction) break;
                        case TWIST: CALL_GENERIC_CODE(phat::bit_tree_pivot_column, phat::twist_reduction) break;
                        case ROW: CALL_GENERIC_CODE(phat::bit_tree_pivot_column, phat::row_reduction) break;
                        case CHUNK: CALL_GENERIC_CODE(phat::bit_tree_pivot_column, phat::chunk_reduction) break;
                        } break;

    case SPARSE_PIVOT:  switch( reduction ) {
                        case STANDARD: CALL_GENERIC_CODE(phat::sparse_pivot_column, phat::standard_reduction) break;
                        case TWIST: CALL_GENERIC_CODE(phat::sparse_pivot_column, phat::twist_reduction) break;
                        case ROW: CALL_GENERIC_CODE(phat::sparse_pivot_column, phat::row_reduction) break;
                        case CHUNK: CALL_GENERIC_CODE(phat::sparse_pivot_column, phat::chunk_reduction) break;
                        } break;
    }
}
