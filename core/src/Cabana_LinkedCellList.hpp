/****************************************************************************
 * Copyright (c) 2018 by the Cabana authors                                 *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the Cabana library. Cabana is distributed under a   *
 * BSD 3-clause license. For the licensing terms see the LICENSE file in    *
 * the top-level directory.                                                 *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#ifndef CABANA_LINKEDCELLLIST_HPP
#define CABANA_LINKEDCELLLIST_HPP

#include <Cabana_Sort.hpp>
#include <Cabana_Slice.hpp>
#include <Cabana_Macros.hpp>
#include <impl/Cabana_CartesianGrid.hpp>

#include <Kokkos_Core.hpp>

namespace Cabana
{
//---------------------------------------------------------------------------//
/*!
  \class LinkedCellList
  \brief Data describing the bin sizes and offsets resulting from a binning
  operation on a 3d regular Cartesian grid.
*/
template<class MemorySpace>
class LinkedCellList
{
  public:

    using memory_space = MemorySpace;
    using KokkosMemorySpace = typename memory_space::kokkos_memory_space;
    using size_type = typename KokkosMemorySpace::size_type;
    using OffsetView = Kokkos::View<size_type*,KokkosMemorySpace>;

    /*!
      \brief Default constructor.
    */
    LinkedCellList()
        : _begin(0)
        , _end(0)
    {}

    /*!
      \brief Slice constructor

      \tparam SliceType Slice type for positions.

      \param positions Slice of positions.

      \param grid_delta Grid sizes in each cardinal direction.

      \param grid_min Grid minimum value in each direction.

      \param grid_max Grid maximum value in each direction.
    */
    template<class SliceType>
    LinkedCellList(
        SliceType positions,
        const typename SliceType::value_type grid_delta[3],
        const typename SliceType::value_type grid_min[3],
        const typename SliceType::value_type grid_max[3],
        typename std::enable_if<(is_slice<SliceType>::value),int>::type * = 0 )
        : _grid( grid_min[0], grid_min[1], grid_min[2],
                 grid_max[0], grid_max[1], grid_max[2],
                 grid_delta[0], grid_delta[1], grid_delta[2] )
        , _begin( 0 )
        , _end( positions.size() )
    {
        build( positions,  grid_delta, grid_min, grid_max );
    }

    /*!
      \brief Slice range constructor

      \tparam SliceType Slice type for positions.

      \param positions Slice of positions.

      \param begin The beginning index of the AoSoA range to sort.

      \param end The end index of the AoSoA range to sort.

      \param grid_delta Grid sizes in each cardinal direction.

      \param grid_min Grid minimum value in each direction.

      \param grid_max Grid maximum value in each direction.
    */
    template<class SliceType>
    LinkedCellList(
        SliceType positions,
        const std::size_t begin,
        const std::size_t end,
        const typename SliceType::value_type grid_delta[3],
        const typename SliceType::value_type grid_min[3],
        const typename SliceType::value_type grid_max[3],
        typename std::enable_if<(is_slice<SliceType>::value),int>::type * = 0 )
        : _grid( grid_min, grid_max, grid_delta ),
          _begin( begin )
        , _end( end )
    {
        build( positions, grid_delta, grid_min, grid_max );
    }

    /*!
      \brief Get the total number of bins.
      \return the total number of bins.
    */
    CABANA_INLINE_FUNCTION
    int totalBins() const
    { return _grid.totalNumCells(); }

    /*!
      \brief Get the number of bins in a given dimension.
      \param dim The dimension to get the number of bins for.
      \return The number of bins.
    */
    CABANA_INLINE_FUNCTION
    int numBin( const int dim ) const
    { return _grid.numBin(dim); }

    /*!
      \brief Given the ijk index of a bin get its cardinal index.
      \param i The i bin index (x).
      \param j The j bin index (y).
      \param k The k bin index (z).
      \return The cardinal bin index.

      Note that the Kokkos sort orders the bins such that the i index moves
      the slowest and the k index mvoes the fastest.
    */
    CABANA_INLINE_FUNCTION
    size_type cardinalBinIndex( const int i, const int j, const int k ) const
    { return _grid.cardinalCellIndex(i,j,k); }

    /*!
      \brief Given the cardinal index of a bin get its ijk indices.
      \param cardinal The cardinal bin index.
      \param i The i bin index (x).
      \param j The j bin index (y).
      \param k The k bin index (z).

      Note that the Kokkos sort orders the bins such that the i index moves
      the slowest and the k index mvoes the fastest.
    */
    CABANA_INLINE_FUNCTION
    void ijkBinIndex( const int cardinal, int& i, int& j, int& k ) const
    {
        _grid.ijkBinIndex( cardinal, i, j, k );
    }

    /*!
      \brief Given a bin get the number of particles it contains.
      \param i The i bin index (x).
      \param j The j bin index (y).
      \param k The k bin index (z).
      \return The number of particles in the bin.
    */
    CABANA_INLINE_FUNCTION
    int binSize( const int i, const int j, const int k ) const
    { return _bin_data.binSize(cardinalBinIndex(i,j,k)); }

    /*!
      \brief Given a bin get the particle index at which it sorts.
      \param i The i bin index (x).
      \param j The j bin index (y).
      \param k The k bin index (z).
      \return The starting particle index of the bin.
    */
    CABANA_INLINE_FUNCTION
    size_type binOffset( const int i, const int j, const int k ) const
    { return _bin_data.binOffset(cardinalBinIndex(i,j,k)); }

    /*!
      \brief Given a local particle id in the binned layout, get the id of the
      particle in the old (unbinned) layout.
      \param particle_id The id of the particle in the binned layout.
      \return The particle id in the old (unbinned) layout.
    */
    CABANA_INLINE_FUNCTION
    size_type permutation( const int particle_id ) const
    { return _bin_data.permutation(particle_id); }

    /*!
      \brief The beginning particle index binned by the linked cell list.
    */
    CABANA_INLINE_FUNCTION
    std::size_t rangeBegin() const
    { return _begin; }

    /*!
      \brief The ending particle index binned by the linked cell list.
    */
    CABANA_INLINE_FUNCTION
    std::size_t rangeEnd() const
    { return _end; }

    /*!
      \brief Get the 1d bin data.
      \return The 1d bin data.
    */
    BinningData<MemorySpace> binningData() const
    { return _bin_data; }

  private:

    template<class SliceType>
    void build( SliceType positions,
                const typename SliceType::value_type grid_delta[3],
                const typename SliceType::value_type grid_min[3],
                const typename SliceType::value_type grid_max[3] )
    {
        // Allocate the binning data. Note that the permutation vector spans
        // only the length of begin-end;
        std::size_t ncell = totalBins();
        Kokkos::View<int*,KokkosMemorySpace> counts( "counts", ncell );
        OffsetView offsets( "offsets", ncell );
        OffsetView permute( "permute", _end - _begin );

        // Count.
        Kokkos::RangePolicy<typename OffsetView::execution_space>
            particle_range( _begin, _end );
        Kokkos::deep_copy( counts, 0 );
        auto cell_count =
            KOKKOS_LAMBDA( const std::size_t p )
            {
                int i, j, k;
                _grid.locatePoint(
                    positions(p,0), positions(p,1), positions(p,2), i , j, k );
                Kokkos::atomic_increment(
                    &counts(_grid.cardinalCellIndex(i,j,k)) );
            };
        Kokkos::parallel_for( "Cabana::LinkedCellList::build::cell_count",
                              particle_range,
                              cell_count );
        Kokkos::fence();

        // Compute offsets.
        Kokkos::RangePolicy<typename OffsetView::execution_space>
            cell_range( 0, ncell );
        auto offset_scan =
            KOKKOS_LAMBDA( const std::size_t c, int& update, const bool final_pass )
            {
                if ( final_pass ) offsets( c ) = update;
                update += counts( c );
            };
        Kokkos::parallel_scan( "Cabana::LinkedCellList::build::offset_scan",
                               cell_range,
                               offset_scan );
        Kokkos::fence();

        // Reset counts.
        Kokkos::deep_copy( counts, 0 );

        // Compute the permutation vector.
        auto create_permute =
            KOKKOS_LAMBDA( const std::size_t p )
            {
                int i, j, k;
                _grid.locatePoint(
                    positions(p,0), positions(p,1), positions(p,2), i , j, k );
                auto cell_id = _grid.cardinalCellIndex(i,j,k);
                int c = Kokkos::atomic_fetch_add( &counts(cell_id), 1 );
                permute( offsets(cell_id) + c ) = p;
            };
        Kokkos::parallel_for( "Cabana::LinkedCellList::build::create_permute",
                              particle_range,
                              create_permute );
        Kokkos::fence();

        // Create the binning data.
        _bin_data = BinningData<MemorySpace>( counts, offsets, permute );
    }

  private:

    BinningData<MemorySpace> _bin_data;
    Impl::CartesianGrid _grid;
    std::size_t _begin;
    std::size_t _end;
};

//---------------------------------------------------------------------------//
// Static type checker.
template<typename >
struct is_linked_cell_list : public std::false_type {};

template<typename MemorySpace>
struct is_linked_cell_list<LinkedCellList<MemorySpace> >
    : public std::true_type {};

template<typename MemorySpace>
struct is_linked_cell_list<const LinkedCellList<MemorySpace> >
    : public std::true_type {};

//---------------------------------------------------------------------------//
/*!
  \brief Given a linked cell list permute an AoSoA.

  \tparam LinkedCellListType The linked cell list type.

  \tparm AoSoA_t The AoSoA type.

  \param permute_vector The permutation vector.

  \param aosoa The AoSoA to permute.
 */
template<class LinkedCellListType, class AoSoA_t>
void permute(
    const LinkedCellListType& linked_cell_list,
    AoSoA_t& aosoa,
    typename std::enable_if<(is_linked_cell_list<LinkedCellListType>::value &&
                             is_aosoa<AoSoA_t>::value),
    int>::type * = 0 )
{
    permute( linked_cell_list.binningData(),
             linked_cell_list.rangeBegin(),
             linked_cell_list.rangeEnd(),
             aosoa );
}

//---------------------------------------------------------------------------//

} // end namespace Cabana

#endif // end CABANA_SORT_HPP