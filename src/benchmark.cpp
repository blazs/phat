/*  Copyright 2013 IST Austria
    Contributed by: Jan Reininghaus

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

//#include <phat/random_access/compute_persistence_pairs.h>

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
#include <phat/random_access/reducers/straight_twist.h>

#include <phat/stack_access/boundary_matrix.h>
#include <phat/stack_access/representations/bit_tree_pivot.h>
#include <phat/stack_access/representations/sparse_pivot.h>
#include <phat/stack_access/representations/full_pivot.h>
#include <phat/stack_access/reducers/standard.h>
#include <phat/stack_access/reducers/twist.h>
#include <phat/stack_access/reducers/straight_twist.h>

#include <phat/auto_reducing/boundary_matrix.h>
#include <phat/auto_reducing/representations/bit_tree_pivot.h>
#include <phat/auto_reducing/representations/sparse_pivot.h>
#include <phat/auto_reducing/representations/full_pivot.h>
#include <phat/auto_reducing/reducers/straight_twist.h>

#include <phat/common/dualize.h>

#include <iostream>
#include <iomanip>

void print_help() {
    std::cerr << "Usage: " << "benchmark " << "[options] input_filename_0 input_filename_1 ... input_filename_N" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << std::endl;
    std::cerr << "--ascii   --  use ascii file format" << std::endl;
    std::cerr << "--binary  --  use binary file format (default)" << std::endl;
    std::cerr << "--help    --  prints this screen" << std::endl;
    std::cerr << "--stack_access   --  use only stack access based algorithms / data structures" << std::endl;
    std::cerr << "--random_access   --  use only random access based algorithms / data structures" << std::endl;
    std::cerr << "--auto_reducing   --  use only auto reducing based algorithms / data structures" << std::endl;
    std::cerr << "--dualize   --  use only dualization approach" << std::endl;
    std::cerr << "--primal   --  use only primal approach" << std::endl;
    std::cerr << "--vector_vector, --vector_set, --vector_list, --full_pivot, --sparse_pivot, --bit_tree_pivot  --  use only a subset of representation data structures for boundary matrices" << std::endl;
    std::cerr << "--standard, --twist, --chunk, --chunk_sequential, --row  --  use only a subset of reduction algorithms" << std::endl;
}

void print_help_and_exit() {
    print_help();
    exit( EXIT_FAILURE );
}

enum Representation_type  {VECTOR_VECTOR, VECTOR_SET, SPARSE_PIVOT, FULL_PIVOT, BIT_TREE_PIVOT, VECTOR_LIST};
enum Algorithm_type  {STANDARD, TWIST, ROW, CHUNK, CHUNK_SEQUENTIAL, STRAIGHT_TWIST};
enum Ansatz_type  {PRIMAL, DUAL};
enum Package_type  {RANDOM_ACCESS, STACK_ACCESS, AUTO_REDUCING};

void parse_command_line( int argc, char** argv, bool& use_binary, std::vector< Representation_type >& representations, std::vector< Algorithm_type >& algorithms
                       , std::vector< Ansatz_type >& ansaetze, std::vector< Package_type >& packages, std::vector< std::string >& input_filenames ) {

    if( argc < 2 ) print_help_and_exit();

    int number_of_options = 0;
    for( int idx = 1; idx < argc; idx++ ) {
        const std::string argument = argv[ idx ];
        if( argument.size() > 2 && argument[ 0 ] == '-' && argument[ 1 ] == '-' ) {
            if( argument == "--ascii" ) use_binary = false;
            else if( argument == "--binary" ) use_binary = true;
            else if( argument == "--vector_vector" ) representations.push_back( VECTOR_VECTOR );
            else if( argument == "--vector_set" ) representations.push_back( VECTOR_SET );
            else if( argument == "--vector_list" ) representations.push_back( VECTOR_LIST );
            else if( argument == "--full_pivot" )  representations.push_back( FULL_PIVOT );
            else if( argument == "--bit_tree_pivot" )  representations.push_back( BIT_TREE_PIVOT );
            else if( argument == "--sparse_pivot" ) representations.push_back( SPARSE_PIVOT );
            else if( argument == "--standard" ) algorithms.push_back( STANDARD );
            else if( argument == "--twist" ) algorithms.push_back( TWIST );
            else if( argument == "--straight_twist" ) algorithms.push_back( STRAIGHT_TWIST );
            else if( argument == "--row" ) algorithms.push_back( ROW );
            else if( argument == "--chunk_sequential" ) algorithms.push_back( CHUNK_SEQUENTIAL );
            else if( argument == "--chunk" ) algorithms.push_back( CHUNK );
            else if( argument == "--primal" ) ansaetze.push_back( PRIMAL );
            else if( argument == "--dual" ) ansaetze.push_back( DUAL );
            else if( argument == "--random_access" ) packages.push_back( RANDOM_ACCESS );
            else if( argument == "--stack_access" ) packages.push_back( STACK_ACCESS );
            else if( argument == "--help" ) print_help_and_exit();
            else print_help_and_exit();
        } else {
            input_filenames.push_back( argument );
        }
    }

    if( representations.empty() == true ) {
        representations.push_back( VECTOR_VECTOR );
        representations.push_back( VECTOR_SET );
        representations.push_back( VECTOR_LIST );
        representations.push_back( FULL_PIVOT );
        representations.push_back( BIT_TREE_PIVOT );
        representations.push_back( SPARSE_PIVOT );
    }

    if( algorithms.empty() == true ) {
        algorithms.push_back( STANDARD );
        algorithms.push_back( TWIST );
        algorithms.push_back( ROW );
        algorithms.push_back( CHUNK );
        algorithms.push_back( CHUNK_SEQUENTIAL );
        algorithms.push_back( STRAIGHT_TWIST );
    }
    
    if( ansaetze.empty() == true ) {
        ansaetze.push_back( PRIMAL );
        ansaetze.push_back( DUAL );
    }

    if( packages.empty() == true ) {
        packages.push_back( RANDOM_ACCESS );
        packages.push_back( STACK_ACCESS );
        packages.push_back( AUTO_REDUCING );
    }
}

template<typename ReducedMatrix, typename Algorithm>
void benchmark( std::string input_filename, bool use_binary, Ansatz_type ansatz ) {

    phat::stack_access::boundary_matrix< phat::stack_access::representations::bit_tree_pivot > matrix;
    bool read_successful = use_binary ? matrix.load_binary( input_filename ) : matrix.load_ascii( input_filename );
   
    if( !read_successful ) {
        std::cerr << std::endl << " Error opening file " << input_filename << std::endl;
        print_help_and_exit();
    }

    ReducedMatrix reduced_matrix;
    Algorithm reduction_algorithm;

    double reduction_timer = -1; 
    if( ansatz == PRIMAL ) {
        std::cout << " primal,";
        reduction_timer = omp_get_wtime();
        reduced_matrix.init( matrix.get_num_cols() );
        reduction_algorithm( matrix, reduced_matrix );
    } else {
        std::cout << " dual,";
        double dualization_timer = omp_get_wtime();
        phat::common::dualize( matrix );
        double dualization_time = omp_get_wtime() - dualization_timer;
        double dualization_time_rounded = floor( dualization_time * 10.0 + 0.5 ) / 10.0;
        std::cout << " Dualization time: " << setiosflags( std::ios::fixed ) << setiosflags( std::ios::showpoint ) << std::setprecision( 1 ) << dualization_time_rounded <<"s,";
        reduction_timer = omp_get_wtime();
        reduced_matrix.init( matrix.get_num_cols() );
        reduction_algorithm( matrix, reduced_matrix );
    }

    double running_time = omp_get_wtime() - reduction_timer;
    double running_time_rounded = floor( running_time * 10.0 + 0.5 ) / 10.0;
    std::cout << " Reduction time: " << setiosflags( std::ios::fixed ) << setiosflags( std::ios::showpoint ) << std::setprecision( 1 ) << running_time_rounded <<"s" << std::endl;
}

#define COMPUTE_RANDOM_ACCESS(Representation) \
    typedef phat::random_access::boundary_matrix< phat::random_access::representations::Representation > ReducedMatrix##Representation;\
    switch( algorithm ) { \
    case STANDARD:         std::cout << input_filename << ", random_access, " << #Representation << ", standard,"; \
                           benchmark< ReducedMatrix##Representation, phat::random_access::reducers::standard >( input_filename, use_binary, ansatz ); \
                           break; \
    case TWIST:            std::cout << input_filename << ", random_access, " << #Representation << ", twist,"; \
                           benchmark< ReducedMatrix##Representation, phat::random_access::reducers::twist >( input_filename, use_binary, ansatz ); \
                           break; \
    case ROW:              std::cout << input_filename << ", random_access, " << #Representation << ", row,"; \
                           benchmark< ReducedMatrix##Representation, phat::random_access::reducers::row >( input_filename, use_binary, ansatz ); \
                           break; \
    case CHUNK:            std::cout << input_filename << ", random_access, " << #Representation << ", chunk,"; \
                           benchmark< ReducedMatrix##Representation, phat::random_access::reducers::chunk >( input_filename, use_binary, ansatz ); \
                           break; \
    case STRAIGHT_TWIST:   std::cout << input_filename << ", random_access, " << #Representation << ", straight_twist,"; \
                           benchmark< ReducedMatrix##Representation, phat::random_access::reducers::straight_twist >( input_filename, use_binary, ansatz ); \
                           break; \
    case CHUNK_SEQUENTIAL: std::cout << input_filename << ", random_access, " << #Representation << ", chunk_sequential,"; \
                           int num_threads = omp_get_max_threads(); \
                           omp_set_num_threads( 1 ); \
                           benchmark< ReducedMatrix##Representation, phat::random_access::reducers::chunk >( input_filename, use_binary, ansatz ); \
                           omp_set_num_threads( num_threads ); \
                           break; \
    };

#define COMPUTE_STACK_ACCESS(Representation) \
    typedef phat::stack_access::boundary_matrix< phat::stack_access::representations::Representation > ReducedMatrix##Representation;\
    switch( algorithm ) { \
    case STANDARD:         std::cout << input_filename << ", stack_access, " << #Representation << ", standard,"; \
                           typedef phat::stack_access::boundary_matrix< phat::stack_access::representations::Representation > ReducedMatrix;\
                           benchmark< ReducedMatrix##Representation, phat::stack_access::reducers::standard >( input_filename, use_binary, ansatz ); \
                           break; \
    case TWIST:            std::cout << input_filename << ", stack_access, " << #Representation << ", twist,"; \
                           typedef phat::stack_access::boundary_matrix< phat::stack_access::representations::Representation > ReducedMatrix;\
                           benchmark< ReducedMatrix##Representation, phat::stack_access::reducers::twist >( input_filename, use_binary, ansatz ); \
                           break; \
    case STRAIGHT_TWIST:   std::cout << input_filename << ", stack_access, " << #Representation << ", straight_twist,"; \
                           benchmark< ReducedMatrix##Representation, phat::stack_access::reducers::straight_twist >( input_filename, use_binary, ansatz ); \
                           break; \
    };

#define COMPUTE_AUTO_REDUCING(Representation) \
    typedef phat::auto_reducing::boundary_matrix< phat::auto_reducing::representations::Representation > ReducedMatrix##Representation;\
    switch( algorithm ) { \
    case STRAIGHT_TWIST:   std::cout << input_filename << ", auto_reducing, " << #Representation << ", straight_twist,"; \
                           benchmark< ReducedMatrix##Representation, phat::auto_reducing::reducers::straight_twist >( input_filename, use_binary, ansatz ); \
                           break; \
    };

int main( int argc, char** argv )
{
    bool use_binary = true; // interpret inputs as binary or ascii files
    std::vector< std::string > input_filenames; // name of file that contains the boundary matrix

    std::vector< Representation_type > representations; // representation class
    std::vector< Algorithm_type > algorithms; // reduction algorithm
    std::vector< Ansatz_type > ansaetze; // primal / dual
    std::vector< Package_type > packages; // stack / random

    parse_command_line( argc, argv, use_binary, representations, algorithms, ansaetze, packages, input_filenames );

    for( int idx_input = 0; idx_input < input_filenames.size(); idx_input++ ) {
        std::string input_filename = input_filenames[ idx_input ];
        for( int idx_algorithm = 0; idx_algorithm < algorithms.size(); idx_algorithm++ ) {
            Algorithm_type algorithm = algorithms[ idx_algorithm ];
            for( int idx_representation = 0; idx_representation < representations.size(); idx_representation++ ) {
                Representation_type cur_representation = representations[ idx_representation ];
                for( int idx_package = 0; idx_package < packages.size(); idx_package++ ) {
                    Package_type cur_package = packages[ idx_package ];
                    for( int idx_ansatz = 0; idx_ansatz < ansaetze.size(); idx_ansatz++ ) {
                        Ansatz_type ansatz = ansaetze[ idx_ansatz ];
                        switch( cur_package ) {
                        case RANDOM_ACCESS: switch( cur_representation ) {
                                            case VECTOR_VECTOR:  COMPUTE_RANDOM_ACCESS(vector_vector) break;
                                            case VECTOR_SET:     COMPUTE_RANDOM_ACCESS(vector_set) break;
                                            case VECTOR_LIST:    COMPUTE_RANDOM_ACCESS(vector_list) break;
                                            case FULL_PIVOT:     COMPUTE_RANDOM_ACCESS(full_pivot) break;
                                            case BIT_TREE_PIVOT: COMPUTE_RANDOM_ACCESS(bit_tree_pivot) break;
                                            case SPARSE_PIVOT:   COMPUTE_RANDOM_ACCESS(sparse_pivot) break;
                                            } break;
                        case STACK_ACCESS:  switch( cur_representation ) {
                                            case BIT_TREE_PIVOT: COMPUTE_STACK_ACCESS(bit_tree_pivot) break;
                                            case FULL_PIVOT:     COMPUTE_STACK_ACCESS(full_pivot) break;
                                            case SPARSE_PIVOT:   COMPUTE_STACK_ACCESS(sparse_pivot) break;
                                            } break;
                        case AUTO_REDUCING: switch( cur_representation ) {
                                            case BIT_TREE_PIVOT: COMPUTE_AUTO_REDUCING(bit_tree_pivot) break;
                                            case FULL_PIVOT:     COMPUTE_AUTO_REDUCING(full_pivot) break;
                                            case SPARSE_PIVOT:   COMPUTE_AUTO_REDUCING(sparse_pivot) break;
                                            } break;
                        }
                    }
                }
            }
        }
    }
}
