//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//



template <typename Iter>
std::set<typename std::iterator_traits<Iter>::value_type>
getDuplicatesInSortedInput(
    Iter begin1, const Iter end1,
    Iter begin2, const Iter end2)
{
    typedef typename std::iterator_traits<Iter>::value_type vtype;

    std::set<vtype> dups;
    while ((begin1 != end1) and (begin2 != end2))
    {
        /// TODO: add an assertion that inputs are sorted?
        if (*begin1 < *begin2)
        {
            begin1++;
        }
        else if (*begin2 < *begin1)
        {
            begin2++;
        }
        else
        {
            dups.insert(*begin1);
            begin1++;
            begin2++;
        }
    }
    return std::move(dups);
}
