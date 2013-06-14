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

#include <phat/random_access/boundary_matrix.h>
#include <phat/random_access/compute_persistence_pairs.h>
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
#include <phat/stack_access/compute_persistence_pairs.h>
#include <phat/stack_access/representations/sparse_pivot.h>
#include <phat/stack_access/representations/full_pivot.h>
#include <phat/stack_access/representations/bit_tree_pivot.h>
#include <phat/stack_access/reducers/standard.h>
#include <phat/stack_access/reducers/twist.h>
#include <phat/stack_access/reducers/straight_twist.h>

#include <phat/auto_reducing/compute_persistence_pairs.h>
#include <phat/auto_reducing/boundary_matrix.h>
#include <phat/auto_reducing/representations/bit_tree_pivot.h>
#include <phat/auto_reducing/representations/sparse_pivot.h>
#include <phat/auto_reducing/representations/full_pivot.h>
#include <phat/auto_reducing/reducers/straight_twist.h>

#define COMPUTE_PAIRS(Package,Representation,Reducer) \
    std::cout << "Running " << #Package << " - " << #Representation << " - " #Reducer " ..." << std::endl;\
    phat::common::persistence_pairs Package##_##Representation##_##Reducer##_##pairs;\
    Package::boundary_matrix< Package::representations::Representation > Package##_##Representation##_##Reducer##_##boundary_matrix = boundary_matrix;\
    Package::compute_persistence_pairs< Package::reducers::Reducer >( Package##_##Representation##_##Reducer##_##pairs, Package##_##Representation##_##Reducer##_##boundary_matrix );

#define COMPARE_PAIRS(FirstPackage,FirstRepresentation,FirstReducer,SecondPackage,SecondRepresentation,SecondReducer) \
    if( FirstPackage##_##FirstRepresentation##_##FirstReducer##_##pairs != SecondPackage##_##SecondRepresentation##_##SecondReducer##_##pairs ) {\
        std::cerr << "Error: " << #FirstPackage << " - " << #FirstRepresentation << " - " #FirstReducer << " and "  << #SecondPackage << " - " << #SecondRepresentation << " - " #SecondReducer << " differ!" << std::endl;\
        error = true;\
    }

int main( int argc, char** argv )
{
    using namespace phat;
    std::string test_data = argc > 1 ? argv[ 1 ] : "torus.bin";

    std::cout << "Reading test data " << test_data << " in binary format ..." << std::endl;
    random_access::boundary_matrix< random_access::representations::full_pivot > boundary_matrix;
    if( !boundary_matrix.load_binary( test_data ) ) {
        std::cerr << "Error: test data " << test_data << " not found!" << std::endl;
        return EXIT_FAILURE;
    }

    bool error = false;
    std::cout << "Comparing representations using (straight-)twist algorithm ..." << std::endl;
    {
        COMPUTE_PAIRS(random_access, sparse_pivot, twist)
        COMPUTE_PAIRS(random_access, full_pivot, twist)
        COMPUTE_PAIRS(random_access, bit_tree_pivot, twist)
        COMPUTE_PAIRS(random_access, vector_vector, twist)
        COMPUTE_PAIRS(random_access, vector_set, twist)
        COMPUTE_PAIRS(random_access, vector_list, twist)
        COMPUTE_PAIRS(stack_access, sparse_pivot, twist)
        COMPUTE_PAIRS(stack_access, full_pivot, twist)
        COMPUTE_PAIRS(stack_access, bit_tree_pivot, twist)
        COMPUTE_PAIRS(auto_reducing, bit_tree_pivot, straight_twist)

        COMPARE_PAIRS(random_access, sparse_pivot, twist, random_access, full_pivot, twist)
        COMPARE_PAIRS(random_access, full_pivot, twist, random_access, bit_tree_pivot, twist)
        COMPARE_PAIRS(random_access, bit_tree_pivot, twist, random_access, vector_vector, twist)
        COMPARE_PAIRS(random_access, vector_vector, twist, random_access, vector_set, twist)
        COMPARE_PAIRS(random_access, vector_set, twist, random_access, vector_list, twist)
        COMPARE_PAIRS(random_access, vector_list, twist, stack_access, sparse_pivot, twist)
        COMPARE_PAIRS(stack_access, sparse_pivot, twist, stack_access, full_pivot, twist)
        COMPARE_PAIRS(stack_access, full_pivot, twist, stack_access, bit_tree_pivot, twist)
        COMPARE_PAIRS(stack_access, bit_tree_pivot, twist, auto_reducing, bit_tree_pivot, straight_twist)
        COMPARE_PAIRS(auto_reducing, bit_tree_pivot, straight_twist, random_access, sparse_pivot, twist)

        if( error ) return EXIT_FAILURE;
        else std::cout << "All results are identical (as they should be)" << std::endl;
    }

    std::cout << "Comparing algorithms using bit_tree_pivot representation ..." << std::endl;
    {
        COMPUTE_PAIRS(random_access, bit_tree_pivot, twist)
        COMPUTE_PAIRS(random_access, bit_tree_pivot, standard)
        COMPUTE_PAIRS(random_access, bit_tree_pivot, chunk)
        COMPUTE_PAIRS(random_access, bit_tree_pivot, row)
        COMPUTE_PAIRS(random_access, bit_tree_pivot, straight_twist)
        COMPUTE_PAIRS(stack_access, bit_tree_pivot, standard)
        COMPUTE_PAIRS(stack_access, bit_tree_pivot, twist)
        COMPUTE_PAIRS(stack_access, bit_tree_pivot, straight_twist)
        COMPUTE_PAIRS(auto_reducing, bit_tree_pivot, straight_twist)

        COMPARE_PAIRS(random_access, bit_tree_pivot, twist, random_access, bit_tree_pivot, standard)
        COMPARE_PAIRS(random_access, bit_tree_pivot, standard, random_access, bit_tree_pivot, chunk)
        COMPARE_PAIRS(random_access, bit_tree_pivot, chunk, random_access, bit_tree_pivot, row)
        COMPARE_PAIRS(random_access, bit_tree_pivot, row, random_access, bit_tree_pivot, straight_twist)
        COMPARE_PAIRS(random_access, bit_tree_pivot, straight_twist, stack_access, bit_tree_pivot, standard)
        COMPARE_PAIRS(stack_access, bit_tree_pivot, standard, stack_access, bit_tree_pivot, twist)
        COMPARE_PAIRS(stack_access, bit_tree_pivot, twist, stack_access, bit_tree_pivot, straight_twist)
        COMPARE_PAIRS(stack_access, bit_tree_pivot, straight_twist, auto_reducing, bit_tree_pivot, straight_twist)
        COMPARE_PAIRS(auto_reducing, bit_tree_pivot, straight_twist, random_access, bit_tree_pivot, twist)

        if( error ) return EXIT_FAILURE;
        else std::cout << "All results are identical (as they should be)" << std::endl;
    }

    std::cout << "Comparing primal and dual approach using random_access and chunk - full_pivot ..." << std::endl;
    {
        common::persistence_pairs primal_pairs;
        random_access::boundary_matrix< random_access::representations::full_pivot > primal_boundary_matrix = boundary_matrix;
        random_access::compute_persistence_pairs< random_access::reducers::chunk >( primal_pairs, primal_boundary_matrix );

        common::persistence_pairs dual_pairs;
        random_access::boundary_matrix< random_access::representations::full_pivot > dual_boundary_matrix = boundary_matrix;
        random_access::compute_persistence_pairs_dualized< random_access::reducers::chunk >( dual_pairs, dual_boundary_matrix );

        if( primal_pairs != dual_pairs ) {
            std::cerr << "Error: primal and dual differ!" << std::endl;
            error = true;
        }

        if( error ) return EXIT_FAILURE;
        else std::cout << "All results are identical (as they should be)" << std::endl;
    }

    std::cout << "Comparing primal and dual approach using stack_access and twist - sparse_pivot ..." << std::endl;
    {
        common::persistence_pairs primal_pairs;
        stack_access::boundary_matrix< stack_access::representations::sparse_pivot > primal_boundary_matrix = boundary_matrix;
        stack_access::compute_persistence_pairs< stack_access::reducers::twist >( primal_pairs, primal_boundary_matrix );

        common::persistence_pairs dual_pairs;
        stack_access::boundary_matrix< stack_access::representations::sparse_pivot > dual_boundary_matrix = boundary_matrix;
        stack_access::compute_persistence_pairs_dualized< stack_access::reducers::twist >( dual_pairs, dual_boundary_matrix );

        if( primal_pairs != dual_pairs ) {
            std::cerr << "Error: primal and dual differ!" << std::endl;
            error = true;
        }

        if( error ) return EXIT_FAILURE;
        else std::cout << "All results are identical (as they should be)" << std::endl;
    }

    std::cout << "Testing normalization ..." << std::endl;
    {
        common::persistence_pairs classic_pairs;
        stack_access::boundary_matrix< stack_access::representations::bit_tree_pivot > classic_boundary_matrix = boundary_matrix;
        stack_access::compute_persistence_pairs< stack_access::reducers::twist >( classic_pairs, classic_boundary_matrix );

        common::persistence_pairs normalized_pairs;
        stack_access::boundary_matrix< stack_access::representations::bit_tree_pivot > normalized_boundary_matrix = boundary_matrix;
        std::vector< phat::index > normalization_map;
        normalized_boundary_matrix.normalize( normalization_map );
        stack_access::compute_persistence_pairs< stack_access::reducers::straight_twist >( normalized_pairs, normalized_boundary_matrix );

        // now normalize classic pairs to compare:
        common::persistence_pairs normalized_classic_pairs;
        for( index idx = 0; idx < classic_pairs.get_num_pairs(); idx++ ) {
            index normalized_birth = normalization_map[ classic_pairs.get_pair( idx ).first ];
            index normalized_death = normalization_map[ classic_pairs.get_pair( idx ).second ];
            normalized_classic_pairs.append_pair( normalized_birth, normalized_death );
        }


        if( normalized_pairs != normalized_classic_pairs ) {
            std::cerr << "Error: normalization bug" << std::endl;
            error = true;
        }

        if( error ) return EXIT_FAILURE;
        else std::cout << "Test passed!" << std::endl;
    }

    std::cout << "Testing random_access vector<vector> interface ..." << std::endl;
    {
        std::vector< std::vector< int > > vector_vector_matrix;
        std::vector< int > vector_dims;
        boundary_matrix.save_vector_vector( vector_vector_matrix, vector_dims );
        random_access::boundary_matrix< random_access::representations::bit_tree_pivot > vector_vector_boundary_matrix;
        vector_vector_boundary_matrix.load_vector_vector( vector_vector_matrix, vector_dims );

        if( vector_vector_boundary_matrix != boundary_matrix ) {
            std::cerr << "Error: [load|save]_vector_vector bug" << std::endl;
            error = true;
        }

        if( error ) return EXIT_FAILURE;
        else std::cout << "Test passed!" << std::endl;
    }

    std::cout << "Testing stack_access vector<vector> interface ..." << std::endl;
    {
        std::vector< std::vector< int > > vector_vector_matrix;
        std::vector< int > vector_dims;
        boundary_matrix.save_vector_vector( vector_vector_matrix, vector_dims );
        stack_access::boundary_matrix< stack_access::representations::bit_tree_pivot > vector_vector_boundary_matrix;
        vector_vector_boundary_matrix.load_vector_vector( vector_vector_matrix, vector_dims );

        if( vector_vector_boundary_matrix != boundary_matrix ) {
            std::cerr << "Error: [load|save]_vector_vector bug" << std::endl;
            error = true;
        }

        if( error ) return EXIT_FAILURE;
        else std::cout << "Test passed!" << std::endl;
    }

    return EXIT_SUCCESS;
}
