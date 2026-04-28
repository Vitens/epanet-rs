//! Reader for EPANET binary `.out` result files.

use std::fs;
const MAGIC_NUMBER: u32 = 516114521;

pub struct EPANETResults {
    pub node_ids: Vec<String>,
    pub link_ids: Vec<String>,
    pub flows: Vec<Vec<f64>>,
    pub heads: Vec<Vec<f64>>,
    pub demands: Vec<Vec<f64>>,
}

/// Helper function to read a fixed length string from a byte array
fn read_fixed_string(bytes: &[u8]) -> String {
    let end = bytes.iter().position(|b| *b == 0).unwrap_or(bytes.len());
    String::from_utf8_lossy(&bytes[..end])
        .trim_end()
        .to_string()
}

/// Read an EPANET outfile and return the results
/// Outfile format: <http://wateranalytics.org/EPANET/_out_file.html>
pub fn read_outfile(filename: &str) -> EPANETResults {
    let data = fs::read(filename).unwrap();

    let n_nodes = u32::from_le_bytes(data[8..12].try_into().unwrap()) as usize;
    let n_tanks = u32::from_le_bytes(data[12..16].try_into().unwrap()) as usize;
    let n_links = u32::from_le_bytes(data[16..20].try_into().unwrap()) as usize;
    let n_pumps = u32::from_le_bytes(data[20..24].try_into().unwrap()) as usize;

    // read epilog
    let epilog = data.len() - 28; // 4 floats + n_periods + warning + magic
    let n_periods = u32::from_le_bytes(data[epilog + 16..epilog + 20].try_into().unwrap()) as usize;
    let magic_number = u32::from_le_bytes(data[epilog + 24..epilog + 28].try_into().unwrap());
    if magic_number != MAGIC_NUMBER {
        panic!("Invalid magic number for epilog");
    }

    let prolog = (884 + 36 * n_nodes + 52 * n_links + 8 * n_tanks) / 4;
    let energy = (28 * n_pumps + 4) / 4;
    let period = (16 * n_nodes + 32 * n_links) / 4;

    let mut flows = vec![vec![0.0; n_links]; n_periods];
    let mut heads = vec![vec![0.0; n_nodes]; n_periods];
    let mut demands = vec![vec![0.0; n_nodes]; n_periods];
    let mut node_ids = Vec::<String>::with_capacity(n_nodes);
    let mut link_ids = Vec::<String>::with_capacity(n_links);

    let node_ids_start = 884usize;
    let node_ids_end = node_ids_start + 32 * n_nodes;
    for chunk in data[node_ids_start..node_ids_end].chunks_exact(32) {
        node_ids.push(read_fixed_string(chunk));
    }
    let link_ids_start = node_ids_end;
    let link_ids_end = link_ids_start + 32 * n_links;
    for chunk in data[link_ids_start..link_ids_end].chunks_exact(32) {
        link_ids.push(read_fixed_string(chunk));
    }

    for (index, chunk) in data.chunks_exact(4).enumerate() {
        let fvalue = f32::from_le_bytes(chunk.try_into().unwrap());
        if index < prolog + energy || index >= prolog + energy + period * n_periods {
            continue;
        }
        let offset = index - prolog - energy;
        let period_index = offset / period;
        let period_offset = offset % period;
        // read demands (1st block of node results)
        if period_offset < n_nodes {
            demands[period_index][period_offset] = fvalue as f64;
            continue;
        }
        // read heads (2nd block of node results)
        if period_offset >= n_nodes && period_offset < n_nodes * 2 {
            heads[period_index][period_offset - n_nodes] = fvalue as f64;
            continue;
        }
        // read flows (1st block of link results)
        if period_offset >= n_nodes * 4 && period_offset < n_nodes * 4 + n_links {
            flows[period_index][period_offset - n_nodes * 4] = fvalue as f64;
            continue;
        }
    }

    return EPANETResults {
        node_ids,
        link_ids,
        flows,
        heads,
        demands,
    };
}
