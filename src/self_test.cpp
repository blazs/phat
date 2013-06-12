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
#include <phat/stack_access/compute_persistence_pairs.h>

#include <phat/random_access/representations/vector_vector.h>
#include <phat/random_access/representations/vector_set.h>
#include <phat/random_access/representations/vector_list.h>
#include <phat/random_access/representations/sparse_pivot.h>
#include <phat/random_access/representations/full_pivot.h>
#include <phat/random_access/representations/bit_tree_pivot.h>

#include <phat/stack_access/representations/bit_tree_pivot.h>
#include <phat/stack_access/reducers/standard.h>

#include <phat/random_access/reducers/twist.h>
#include <phat/random_access/reducers/standard.h>
#include <phat/random_access/reducers/row.h>
#include <phat/random_access/reducers/chunk.h>

int main( int argc, char** argv )
{
    std::string test_data = argc > 1 ? argv[ 1 ] : "torus.bin";

    typedef phat::random_access::representations::sparse_pivot Sparse;
    typedef phat::random_access::representations::full_pivot Full;
    typedef phat::random_access::representations::bit_tree_pivot BitTree;
    typedef phat::random_access::representations::vector_vector Vec_vec;
    typedef phat::random_access::representations::vector_set Vec_set;
    typedef phat::random_access::representations::vector_list Vec_list;

    std::cout << "Reading test data " << test_data << " in binary format ..." << std::endl;
    phat::random_access::boundary_matrix< Full > boundary_matrix;
    if( !boundary_matrix.load_binary( test_data ) ) {
        std::cerr << "Error: test data " << test_data << " not found!" << std::endl;
        return EXIT_FAILURE;
    }

    bool error = false;
    std::cout << "Comparing representations using Chunk algorithm ..." << std::endl;
    {
        std::cout << "Running Chunk - Sparse ..." << std::endl;
        phat::common::persistence_pairs sparse_pairs;
        phat::random_access::boundary_matrix< Sparse > sparse_boundary_matrix = boundary_matrix;
        phat::random_access::compute_persistence_pairs< phat::random_access::reducers::chunk >( sparse_pairs, sparse_boundary_matrix );

        std::cout << "Running Chunk - Full ..." << std::endl;
        phat::common::persistence_pairs full_pairs;
        phat::random_access::boundary_matrix< Full > full_boundary_matrix = boundary_matrix;
        phat::random_access::compute_persistence_pairs< phat::random_access::reducers::chunk >( full_pairs, full_boundary_matrix );

        std::cout << "Running Chunk - BitTree ..." << std::endl;
        phat::common::persistence_pairs bit_tree_pairs;
        phat::random_access::boundary_matrix< BitTree > bit_tree_boundary_matrix = boundary_matrix;
        phat::random_access::compute_persistence_pairs< phat::random_access::reducers::chunk >( bit_tree_pairs, bit_tree_boundary_matrix );

        std::cout << "Running Chunk - Vec_vec ..." << std::endl;
        phat::common::persistence_pairs vec_vec_pairs;
        phat::random_access::boundary_matrix< Vec_vec > vec_vec_boundary_matrix = boundary_matrix;
        phat::random_access::compute_persistence_pairs< phat::random_access::reducers::chunk >( vec_vec_pairs, vec_vec_boundary_matrix );

        std::cout << "Running Chunk - Vec_set ..." << std::endl;
        phat::common::persistence_pairs vec_set_pairs;
        phat::random_access::boundary_matrix< Vec_set > vec_set_boundary_matrix = boundary_matrix;
        phat::random_access::compute_persistence_pairs< phat::random_access::reducers::chunk >( vec_set_pairs, vec_set_boundary_matrix );

        std::cout << "Running Chunk - Vec_list ..." << std::endl;
        phat::common::persistence_pairs vec_list_pairs;
        phat::random_access::boundary_matrix< Vec_list > vec_list_boundary_matrix = boundary_matrix;
        phat::random_access::compute_persistence_pairs< phat::random_access::reducers::chunk >( vec_list_pairs, vec_list_boundary_matrix );

        if( sparse_pairs != full_pairs ) {
            std::cerr << "Error: sparse and full differ!" << std::endl;
            error = true;
        }
        if( full_pairs != vec_vec_pairs ) {
            std::cerr << "Error: full and vec_vec differ!" << std::endl;
            error = true;
        }
        if( vec_vec_pairs != vec_set_pairs ) {
            std::cerr << "Error: vec_vec and vec_set differ!" << std::endl;
            error = true;
        }
        if( vec_set_pairs != bit_tree_pairs ) {
            std::cerr << "Error: vec_set and bit_tree differ!" << std::endl;
            error = true;
        }
         if( bit_tree_pairs != vec_list_pairs ) {
            std::cerr << "Error: bit_tree and vec_list differ!" << std::endl;
            error = true;
        }
        if( vec_list_pairs != sparse_pairs ) {
            std::cerr << "Error: vec_list and sparse differ!" << std::endl;
            error = true;
        }

        if( error ) return EXIT_FAILURE;
        else std::cout << "All results are identical (as they should be)" << std::endl;
    }

    std::cout << "Comparing algorithms using BitTree representation ..." << std::endl;
    {
        std::cout << "Running Twist - BitTree ..." << std::endl;
        phat::common::persistence_pairs twist_pairs;
        phat::random_access::boundary_matrix< BitTree > twist_boundary_matrix = boundary_matrix;
        phat::random_access::compute_persistence_pairs< phat::random_access::reducers::twist >( twist_pairs, twist_boundary_matrix );

        std::cout << "Running Standard - BitTree ..." << std::endl;
        phat::common::persistence_pairs std_pairs;
        phat::random_access::boundary_matrix< BitTree > std_boundary_matrix = boundary_matrix;
        phat::random_access::compute_persistence_pairs< phat::random_access::reducers::standard >( std_pairs, std_boundary_matrix );

        std::cout << "Running Chunk - BitTree ..." << std::endl;
        phat::common::persistence_pairs chunk_pairs;
        phat::random_access::boundary_matrix< BitTree > chunk_boundary_matrix = boundary_matrix;
        phat::random_access::compute_persistence_pairs< phat::random_access::reducers::chunk >( chunk_pairs, chunk_boundary_matrix );

        std::cout << "Running Row - BitTree ..." << std::endl;
        phat::common::persistence_pairs row_pairs;
        phat::random_access::boundary_matrix< BitTree > row_boundary_matrix = boundary_matrix;
        phat::random_access::compute_persistence_pairs< phat::random_access::reducers::row >( row_pairs, row_boundary_matrix );

        if( twist_pairs != std_pairs ) {
            std::cerr << "Error: twist and standard differ!" << std::endl;
            error = true;
        }
        if( std_pairs != chunk_pairs ) {
            std::cerr << "Error: standard and chunk differ!" << std::endl;
            error = true;
        }
        if( chunk_pairs != row_pairs ) {
            std::cerr << "Error: chunk and row differ!" << std::endl;
            error = true;
        }
        if( row_pairs != twist_pairs ) {
            std::cerr << "Error: row and twist differ!" << std::endl;
            error = true;
        }

        if( error ) return EXIT_FAILURE;
        else std::cout << "All results are identical (as they should be)" << std::endl;
    }

    std::cout << "Comparing random_access to stack_access using Standard - BitTree ..." << std::endl;
    {
        phat::common::persistence_pairs random_access_pairs;
        phat::random_access::boundary_matrix< BitTree > random_access_boundary_matrix = boundary_matrix;
        phat::random_access::compute_persistence_pairs< phat::random_access::reducers::standard >( random_access_pairs, random_access_boundary_matrix );

        phat::common::persistence_pairs stack_access_pairs;
        phat::stack_access::boundary_matrix< phat::stack_access::representations::bit_tree_pivot > stack_access_boundary_matrix = boundary_matrix;
        phat::stack_access::compute_persistence_pairs< phat::stack_access::reducers::standard >( stack_access_pairs, stack_access_boundary_matrix );

        if( random_access_pairs != stack_access_pairs ) {
            std::cerr << "Error: random_access and stack_access differ!" << std::endl;
            error = true;
        }

        if( error ) return EXIT_FAILURE;
        else std::cout << "All results are identical (as they should be)" << std::endl;
    }

    std::cout << "Comparing primal and dual approach using Chunk - Full ..." << std::endl;
    {
        phat::common::persistence_pairs primal_pairs;
        phat::random_access::boundary_matrix< Full > primal_boundary_matrix = boundary_matrix;
        phat::random_access::compute_persistence_pairs< phat::random_access::reducers::chunk >( primal_pairs, primal_boundary_matrix );

        phat::common::persistence_pairs dual_pairs;
        phat::random_access::boundary_matrix< Full > dual_boundary_matrix = boundary_matrix;
        phat::random_access::compute_persistence_pairs_dualized< phat::random_access::reducers::chunk >( dual_pairs, dual_boundary_matrix );

        if( primal_pairs != dual_pairs ) {
            std::cerr << "Error: primal and dual differ!" << std::endl;
            error = true;
        }

        if( error ) return EXIT_FAILURE;
        else std::cout << "All results are identical (as they should be)" << std::endl;
    }

    std::cout << "Testing vector<vector> interface ..." << std::endl;
    {
        std::vector< std::vector< int > > vector_vector_matrix;
        std::vector< int > vector_dims;
        boundary_matrix.save_vector_vector( vector_vector_matrix, vector_dims );
        phat::random_access::boundary_matrix< BitTree > vector_vector_boundary_matrix;
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
