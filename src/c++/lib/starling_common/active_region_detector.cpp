#include "active_region_detector.hh"
#include <cstdio>

void
active_region_detector::insert_mismatch_position(pos_t pos)
{
//    printf("\tMISMATCH\t%d\n", pos);
    add_count(pos, 1);
}

void
active_region_detector::insert_indel_position(pos_t pos)
{
//    printf("\tINDEL\t%d\n", pos);
    add_count(pos, 4);
}

static int index = 0;

void
active_region_detector::update_start_position(const pos_t pos)
{
//    printf("START\t%d\n", pos);
    for (int i(_buffer_start_pos); i<pos; ++i) reset_count(i);

    _buffer_start_pos = pos;

    if (_active_regions.empty()) return;

    active_region front =  _active_regions.front();
    if (front.get_end() < _buffer_start_pos)
    {
        printf("%d\t%d\t%d\t%d\t%d\n", ++index, front.get_start()+1, front.get_end()+1, front.get_num_variants(), front.get_length());
        _active_regions.pop_front();
    }
}

void
active_region_detector::update_end_position(const pos_t pos)
{
//    printf("END\t%d\n", pos);
    if (!is_candidate_variant(pos)) return;

    if (pos - _active_region_start_pos  < (int)_max_detection_window_size)
    {
        // extend the existing active region to include this position
        ++_num_variants;
    }
    else
    {
        // close existing active region
        if (_num_variants >= _min_num_variants_per_region)
        {
            _active_regions.push_back(active_region(_active_region_start_pos, _prev_variant_pos, _num_variants));
        }
        // start new active region
        _active_region_start_pos = pos;
        _num_variants = 1;
    }
    _prev_variant_pos = pos;

    // TODO: remaining
}

bool
active_region_detector::is_candidate_variant(const pos_t pos) const
{
//    if (get_count(pos) >= _min_num_mismatches_per_position)
//        printf("*** %d, %d\n", (pos+1), get_count(pos));
    return get_count(pos) >= _min_num_mismatches_per_position;
}

std::list<active_region>&
active_region_detector::get_active_regions()
{
    return _active_regions;
}
