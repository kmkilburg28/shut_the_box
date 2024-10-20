use wasm_bindgen::prelude::*;
use std::{collections::HashSet, cmp::min};


fn combinations(n: u64, k: u64) -> f64 {
    if k > n {
        0.0
    } else {
        let mut comb = 1.0;
        for j in 0..min(k, n - k) {
            comb /= (j + 1) as f64;
            comb *= (n - j) as f64;
        }
        comb
    }
}

fn fill_out_number_masks(min_sum: u32, max_sum: u32, game_size: u32, indices_masks: &Vec<u32>) -> Vec<HashSet<u32>> {
    let mut number_maps: Vec<HashSet<u32>> = vec![HashSet::new(); (max_sum+1 - min_sum) as usize];
    let number_maps_below_min: Vec<HashSet<u32>> = (0..(min_sum-1) as usize).map(|i| HashSet::from([indices_masks[i]])).collect();

    for base_sum in min_sum..max_sum+1 {
        // let number_map = &mut number_maps[(base_sum-min_sum) as usize];
        let cur_ind = (base_sum-min_sum) as usize;
        if base_sum <= game_size {
            number_maps[cur_ind].insert(indices_masks[(base_sum-1) as usize]);
        }

        for low_num in 1..(base_sum/2)+1 {
            let high_num = base_sum - low_num;

            let low_num_masks = if low_num < min_sum { &number_maps_below_min[(low_num-1) as usize] } else { &number_maps[(low_num-min_sum) as usize] };
            let high_num_masks = if high_num < min_sum { &number_maps_below_min[(high_num-1) as usize] } else { &number_maps[(high_num-min_sum) as usize] };

            // Duplicate HashSet to save masks into to get around rust's memory safety around the number_maps vector
            let mut new_masks:HashSet<u32> = HashSet::new();
            for low_mask in low_num_masks {
                for high_mask in high_num_masks {
                    let or_mask = low_mask | high_mask;
                    let xor_compare = low_mask ^ high_mask;
                    if or_mask == xor_compare {
                        new_masks.insert(or_mask);
                    }
                }
            }
            number_maps[cur_ind].extend(new_masks);
        }
    }

    number_maps
}

fn get_next_states(min_sum: u32, result_sum: u32, gamestate: u32, number_masks: &Vec<HashSet<u32>>) -> HashSet<u32> {
    let mut next_states: HashSet<u32> = HashSet::new();
    for mask in &number_masks[(result_sum-min_sum) as usize] {
        if gamestate & *mask == *mask {
            next_states.insert(gamestate & !mask);
        }
    }
    next_states
}

fn fill_out_states_map(min_sum: u32, max_sum: u32, num_states: u32, number_masks: &Vec<HashSet<u32>>) -> Vec<Vec<HashSet<u32>>> {
    // let mut states_map:Vec<Vec<HashSet<u32>>> = (0..num_states).map(|_| (0..max_sum-min_sum+1).map(|_| HashSet::new()).collect::<Vec<HashSet<u32>>>()).collect();
    let mut states_map:Vec<Vec<HashSet<u32>>> = vec![vec![HashSet::new(); (max_sum-min_sum+1) as usize]; num_states as usize];
    for state in 0..num_states {
        for sum in min_sum..max_sum+1 {
            states_map[state as usize][(sum-min_sum) as usize] = get_next_states(min_sum, sum, state, number_masks);
        }
    }
    states_map
}

fn calculate_sum_odds(die_count: u32, die_size: u32) -> Vec<f64> {
    let min_sum = die_count;
    let max_sum = die_count*die_size;
    let mut sum_odds:Vec<f64> = vec![0.0; (max_sum-min_sum+1) as usize];

    let permutations_count = u32::pow(die_size, die_count);
    // Can be optimized by halving the calculations and reflecting it over
    // Thanks to https://www.quora.com/Is-there-a-formula-to-calculate-the-probability-of-the-sum-of-x-dice-being-than-y
    for sum in min_sum..max_sum+1 {
        let ind = (sum - min_sum) as usize;
        for i in 0..((sum-die_count)/die_size)+1 {
            sum_odds[ind] += (if i % 2 == 0 {1.0} else {-1.0}) * combinations(die_count as u64,i as u64) * combinations((sum-(die_size*i)-1) as u64, (die_count-1) as u64);
        }
        sum_odds[ind] /= permutations_count as f64;
    }
    sum_odds
}

fn generate_gameboard(game_size: u32) -> u32 {
    u32::pow(2,game_size) - 1
}

fn state_to_num_sum(game_size: u32, gamestate: u32, indices_masks: &Vec<u32>) -> u32 {
    (1..game_size+1).filter(|i| indices_masks[(i-1) as usize] & gamestate != 0).sum()
}

fn score_states(die_count: u32, die_size: u32, game_size: u32, indices_masks: &Vec<u32>, states_map: &Vec<Vec<HashSet<u32>>>) -> Vec<f64> {
    let sum_odds = calculate_sum_odds(die_count, die_size);
    let max_sum_score:f64 = state_to_num_sum(game_size, generate_gameboard(game_size), indices_masks).into();

    let mut state_scores:Vec<f64> = vec![max_sum_score; states_map.len()];
    for (state,state_map) in states_map.iter().enumerate() {
        let mut avg_score: f64 = 0.0;
        // Iterate through state sets resulting from a specific sum
        for (i,sum_set) in state_map.iter().enumerate() {
            let mut sum_set_score:f64 = max_sum_score;
            for to_state in sum_set {
                let to_state_score = state_scores[*to_state as usize];
                if to_state_score < sum_set_score {
                    sum_set_score = to_state_score;
                }
            }
            if sum_set.len() == 0 {
                sum_set_score = state_to_num_sum(game_size, state as u32, indices_masks).into();
            }

            avg_score += sum_odds[i]*sum_set_score;
        }
        state_scores[state] = avg_score;
    }
    state_scores
}

#[wasm_bindgen]
pub fn run(die_count: u32, die_size: u32, game_size: u32) -> Vec<f64> {
    let min_sum = die_count;
    let max_sum = die_count * die_size;
    let num_states = u32::pow(2, game_size);
    let indices_masks:Vec<u32> = (0..game_size).map(|x| u32::pow(2,x)).collect();

    let number_masks = fill_out_number_masks(min_sum, max_sum, game_size, &indices_masks);
    let states_map = fill_out_states_map(min_sum, max_sum, num_states, &number_masks);
    let state_scores = score_states(die_count, die_size, game_size, &indices_masks, &states_map);
    for i in 0..state_scores.len() {
        println!("{0}: {1:b} {2}", i, i, state_scores[i]);
    }
    state_scores
}
