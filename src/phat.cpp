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

#include <phat/random_access/compute_persistence_pairs.h>
#include <phat/random_access/boundary_matrix.h>
#include <phat/random_access/representations/vector_vector.h>
#include <phat/random_access/representations/vector_set.h>
#include <phat/random_access/representations/vector_list.h>
#include <phat/random_access/representations/sparse_pivot.h>
#include <phat/random_access/representations/full_pivot.h>
#include <phat/random_access/representations/bit_tree_pivot.h>
#include <phat/random_access/reducers/twist.h>
#include <phat/random_access/reducers/standard.h>
#include <phat/random_access/reducers/row.h>
#include <phat/random_access/reducers/chunk.h>

#include <phat/stack_access/compute_persistence_pairs.h>
#include <phat/stack_access/boundary_matrix.h>
#include <phat/stack_access/representations/bit_tree_pivot.h>
#include <phat/stack_access/representations/sparse_pivot.h>
#include <phat/stack_access/representations/full_pivot.h>
#include <phat/stack_access/representations/bit_tree_pivot.h>
#include <phat/stack_access/reducers/standard.h>

#include <phat/common/dualize.h>

enum Representation_type  {VECTOR_VECTOR, VECTOR_SET, SPARSE_PIVOT, FULL_PIVOT, BIT_TREE_PIVOT, VECTOR_LIST};
enum Algorithm_type  {STANDARD, TWIST, ROW, CHUNK, CHUNK_SEQUENTIAL };
enum Package_type  {RANDOM_ACCESS, STACK_ACCESS};

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
    std::cerr << "--stack_access, --random_access   --  selects a package (default is '--random_access')" << std::endl;
    std::cerr << "--vector_vector, --vector_set, --vector_list, --full_pivot, --sparse_pivot, --bit_tree_pivot  --  selects a representation data structure for boundary matrices (default is '--bit_tree_pivot')" << std::endl;
    std::cerr << "--standard, --twist, --chunk, --chunk_sequential, --row  --  selects a reduction algorithm (default is '--twist')" << std::endl;
}

void print_help_and_exit() {
    print_help();
    exit( EXIT_FAILURE );
}

void parse_command_line( int argc, char** argv, bool& use_binary, Package_type& package, Representation_type& representation, Algorithm_type& algorithm,
                         std::string& input_filename, std::string& output_filename, bool& verbose, bool& dualize) {

    if( argc < 3 ) print_help_and_exit();

    input_filename = argv[ argc - 2 ];
    output_filename = argv[ argc - 1 ];

    for( int idx = 1; idx < argc - 2; idx++ ) {
        const std::string option = argv[ idx ];

        if( option == "--ascii" ) use_binary = false;
        else if( option == "--binary" ) use_binary = true;
        else if( option == "--dualize" ) dualize = true;
        else if( option == "--random_access" ) package = RANDOM_ACCESS;
        else if( option == "--stack_access" ) package = STACK_ACCESS;
        else if( option == "--vector_vector" ) representation = VECTOR_VECTOR;
        else if( option == "--vector_set" ) representation = VECTOR_SET;
        else if( option == "--vector_list" ) representation = VECTOR_LIST;
        else if( option == "--full_pivot" )  representation = FULL_PIVOT;
        else if( option == "--bit_tree_pivot" )  representation = BIT_TREE_PIVOT;
        else if( option == "--sparse_pivot" ) representation = SPARSE_PIVOT;
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

template< typename BoundaryMatrix >
void load_boundary_matrix( BoundaryMatrix& matrix, std::string input_filename, bool use_binary, bool verbose, bool dualize ) {
    bool read_successful;

    double read_timer = omp_get_wtime();
    if( use_binary ) {
        LOG( "Reading input file " << input_filename << " in binary mode" )
        read_successful = matrix.load_binary( input_filename );
    } else {
        LOG( "Reading input file " << input_filename << " in ascii mode" )
        read_successful = matrix.load_ascii( input_filename );
    }
    double read_time = omp_get_wtime() - read_timer;
    double read_time_rounded = floor( read_time * 10.0 + 0.5 ) / 10.0;
    LOG( "Reading input file took " << setiosflags( std::ios::fixed ) << setiosflags( std::ios::showpoint ) << std::setprecision( 1 ) << read_time_rounded <<"s" )

    if( !read_successful ) {
        std::cerr << "Error opening file " << input_filename << std::endl;
        print_help_and_exit();
    }

    if( dualize ) {
        double dualize_timer = omp_get_wtime();
        LOG( "Dualizing ..." )
        phat::common::dualize ( matrix );
        double dualize_time = omp_get_wtime() - dualize_timer;
        double dualize_time_rounded = floor( dualize_time * 10.0 + 0.5 ) / 10.0;
        LOG( "Dualizing took " << setiosflags( std::ios::fixed ) << setiosflags( std::ios::showpoint ) << std::setprecision( 1 ) << dualize_time_rounded <<"s" )
    }
}

void save_persistence_pairs( phat::common::persistence_pairs& pairs, std::string output_filename, bool use_binary, bool verbose ) {
    double write_timer = omp_get_wtime();
    if( use_binary ) {
        LOG( "Writing output file " << output_filename << " in binary mode ..." )
        pairs.save_binary( output_filename );
    } else {
        LOG( "Writing output file " << output_filename << " in ascii mode ..." )
        pairs.save_ascii( output_filename );
    }
    double write_time = omp_get_wtime() - write_timer;
    double write_time_rounded = floor( write_time * 10.0 + 0.5 ) / 10.0;
    LOG( "Writing output file took " << setiosflags( std::ios::fixed ) << setiosflags( std::ios::showpoint ) << std::setprecision( 1 ) << write_time_rounded <<"s" )

}

template<typename Representation, typename Algorithm>
void compute_pairing_random_access( std::string input_filename, std::string output_filename, bool use_binary, bool verbose, bool dualize ) {
    
    phat::random_access::boundary_matrix< Representation > matrix;
    load_boundary_matrix( matrix, input_filename, use_binary, verbose, dualize );
    
    double pairs_timer = omp_get_wtime();
    phat::common::persistence_pairs pairs;
    LOG( "Computing persistence pairs ..." )
    phat::random_access::compute_persistence_pairs< Algorithm >( pairs, matrix );
    if( dualize ) phat::common::dualize_persistence_pairs( pairs, matrix.get_num_cols() );
    double pairs_time = omp_get_wtime() - pairs_timer;
    double pairs_time_rounded = floor( pairs_time * 10.0 + 0.5 ) / 10.0;
    LOG( "Computing persistence pairs took " << setiosflags( std::ios::fixed ) << setiosflags( std::ios::showpoint ) << std::setprecision( 1 ) << pairs_time_rounded <<"s" )

    save_persistence_pairs( pairs, output_filename, use_binary, verbose );
}

template<typename Representation, typename Algorithm>
void compute_pairing_stack_access( std::string input_filename, std::string output_filename, bool use_binary, bool verbose, bool dualize ) {
    
    phat::stack_access::boundary_matrix< Representation > matrix;
    load_boundary_matrix( matrix, input_filename, use_binary, verbose, dualize );
    
    double pairs_timer = omp_get_wtime();
    phat::common::persistence_pairs pairs;
    LOG( "Computing persistence pairs ..." )
    phat::stack_access::compute_persistence_pairs< Algorithm >( pairs, matrix );
    if( dualize ) phat::common::dualize_persistence_pairs( pairs, matrix.get_num_cols() );
    double pairs_time = omp_get_wtime() - pairs_timer;
    double pairs_time_rounded = floor( pairs_time * 10.0 + 0.5 ) / 10.0;
    LOG( "Computing persistence pairs took " << setiosflags( std::ios::fixed ) << setiosflags( std::ios::showpoint ) << std::setprecision( 1 ) << pairs_time_rounded <<"s" )

    save_persistence_pairs( pairs, output_filename, use_binary, verbose );
}

#define CALL_COMPUTE_PAIRING_RANDOM_ACCESS(Representation) \
    switch( algorithm ) { \
    case STANDARD: compute_pairing_random_access< phat::random_access::representations::Representation, phat::random_access::reducers::standard > ( input_filename, output_filename, use_binary, verbose, dualize ); break; \
    case TWIST: compute_pairing_random_access< phat::random_access::representations::Representation, phat::random_access::reducers::twist > ( input_filename, output_filename, use_binary, verbose, dualize ); break; \
    case ROW: compute_pairing_random_access< phat::random_access::representations::Representation, phat::random_access::reducers::row >( input_filename, output_filename, use_binary, verbose, dualize ); break; \
    case CHUNK: compute_pairing_random_access< phat::random_access::representations::Representation, phat::random_access::reducers::chunk >( input_filename, output_filename, use_binary, verbose, dualize ); break; \
    case CHUNK_SEQUENTIAL: { int num_threads = omp_get_max_threads(); \
                           omp_set_num_threads( 1 ); \
                           compute_pairing_random_access< phat::random_access::representations::Representation, phat::random_access::reducers::chunk >( input_filename, output_filename, use_binary, verbose, dualize ); break; \
                           omp_set_num_threads( num_threads ); \
                           } break; \
    default: std::cerr << "Invalid package-representation-algorithm combination. Please refer to documentation or source code for valid combinations." << std::endl; break;\
    }

#define CALL_COMPUTE_PAIRING_STACK_ACCESS(Representation) \
    switch( algorithm ) { \
    case STANDARD: compute_pairing_stack_access< phat::stack_access::representations::Representation, phat::stack_access::reducers::standard > ( input_filename, output_filename, use_binary, verbose, dualize ); break; \
    default: std::cerr << "Invalid package-representation-algorithm combination. Please refer to documentation or source code for valid combinations." << std::endl;\
    }

int main( int argc, char** argv )
{
    bool use_binary = true; // interpret input as binary or ascii file
    Package_type package = RANDOM_ACCESS; // package to use
    Representation_type representation = BIT_TREE_PIVOT; // representation class
    Algorithm_type algorithm = TWIST; // reduction algorithm
    std::string input_filename; // name of file that contains the boundary matrix
    std::string output_filename; // name of file that will contain the persistence pairs
    bool verbose = false; // print timings / info
    bool dualize = false; // toggle for dualization approach

    parse_command_line( argc, argv, use_binary, package, representation, algorithm, input_filename, output_filename, verbose, dualize );

    switch( package ) {
    case RANDOM_ACCESS: 
        switch( representation ) {
        case VECTOR_VECTOR: CALL_COMPUTE_PAIRING_RANDOM_ACCESS(vector_vector) break;
        case VECTOR_SET: CALL_COMPUTE_PAIRING_RANDOM_ACCESS(vector_set) break;
        case VECTOR_LIST: CALL_COMPUTE_PAIRING_RANDOM_ACCESS(vector_list) break;
        case FULL_PIVOT: CALL_COMPUTE_PAIRING_RANDOM_ACCESS(full_pivot) break;
        case BIT_TREE_PIVOT: CALL_COMPUTE_PAIRING_RANDOM_ACCESS(bit_tree_pivot) break;
        case SPARSE_PIVOT: CALL_COMPUTE_PAIRING_RANDOM_ACCESS(sparse_pivot) break;
        default: std::cerr << "Invalid package-representation-algorithm combination. Please refer to documentation or source code for valid combinations." << std::endl;
        } break;
    case STACK_ACCESS:
        switch( representation ) {
        case FULL_PIVOT: CALL_COMPUTE_PAIRING_STACK_ACCESS(full_pivot) break;
        case BIT_TREE_PIVOT: CALL_COMPUTE_PAIRING_STACK_ACCESS(bit_tree_pivot) break;
        case SPARSE_PIVOT: CALL_COMPUTE_PAIRING_STACK_ACCESS(sparse_pivot) break;
        default: std::cerr << "Invalid package-representation-algorithm combination. Please refer to documentation or source code for valid combinations" << std::endl;
        } break;
    }
}
