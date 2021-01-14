/****************************************************************************
 * Copyright (c) 2018-2020 by the Cabana authors                            *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the Cabana library. Cabana is distributed under a   *
 * BSD 3-clause license. For the licensing terms see the LICENSE file in    *
 * the top-level directory.                                                 *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#include <Cajita_UniformDimPartitioner.hpp>

namespace Cajita
{
//---------------------------------------------------------------------------//
std::array<int, 3>
UniformDimPartitioner::ranksPerDimension( MPI_Comm comm,
                                          const std::array<int, 3>& ) const
{
    int comm_size;
    MPI_Comm_size( comm, &comm_size );

    std::array<int, 3> ranks_per_dim = { 0, 0, 0 };
    MPI_Dims_create( comm_size, 3, ranks_per_dim.data() );

    return ranks_per_dim;
}

//---------------------------------------------------------------------------//
std::array<int,3>
UniformDimPartitioner::ownedCellsPerDimension(
    MPI_Comm cart_comm,
    const std::array<int,3>& ranks_per_dim,
    const std::array<int, 3>& global_cells_per_dim ) const
{
    // Get the Cartesian topology index of this rank.
    std::array<int,3> cart_rank;
    int linear_rank;
    MPI_Comm_rank( cart_comm, &linear_rank );
    MPI_Cart_coords( cart_comm, linear_rank, 3, cart_rank.data() );

    // Get the cells per dimension and the remainder.
    std::array<int, 3> cells_per_dim;
    std::array<int, 3> dim_remainder;
    for ( int d = 0; d < 3; ++d )
    {
        cells_per_dim[d] = global_num_cell[d] / ranks_per_dim[d];
        dim_remainder[d] = global_num_cell[d] % ranks_per_dim[d];
    }

    // Compute the number of local cells in this rank in each dimension.
    std::array<int,3> owned_num_cell
    for ( int d = 0; d < 3; ++d )
    {
        owned_num_cell[d] = cells_per_dim[d];
        if ( dim_remainder[d] > _cart_rank[d] )
            ++owned_num_cell[d];
    }

    return owned_num_cell;
}

//---------------------------------------------------------------------------//

} // end namespace Cajita
