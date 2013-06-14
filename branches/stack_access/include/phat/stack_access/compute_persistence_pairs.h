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

#include <phat/common/persistence_pairs.h>
#include <phat/common/dualize.h>

#include <phat/stack_access/boundary_matrix.h>
#include <phat/stack_access/reducers/standard.h>

namespace phat { namespace stack_access {

    template< typename ReductionAlgorithm, typename Representation >
    void compute_persistence_pairs( common::persistence_pairs& pairs, boundary_matrix< Representation>& boundary_matrix ) {
        ReductionAlgorithm reduce;
        phat::stack_access::boundary_matrix< Representation> reduced_matrix;
        reduced_matrix.init( boundary_matrix.get_num_cols() );
        reduce( boundary_matrix, reduced_matrix );
        pairs.read_off_pairs( reduced_matrix );
    }
    
    template< typename ReductionAlgorithm, typename Representation >
    void compute_persistence_pairs_dualized( common::persistence_pairs& pairs, boundary_matrix< Representation >& boundary_matrix ) {
        common::dualize( boundary_matrix );
        compute_persistence_pairs< ReductionAlgorithm >( pairs, boundary_matrix );
        common::dualize_persistence_pairs( pairs, boundary_matrix.get_num_cols() );
    }
    
    template< typename Representation >
    void compute_persistence_pairs( common::persistence_pairs& pairs, boundary_matrix< Representation >& boundary_matrix ) {
        compute_persistence_pairs< reducers::standard >( pairs, boundary_matrix );
    }
    
    
    template< typename Representation >
    void compute_persistence_pairs_dualized( common::persistence_pairs& pairs, boundary_matrix< Representation >& boundary_matrix ) {
        compute_persistence_pairs_dualized< reducers::standard >( pairs, boundary_matrix );
    }
} }



