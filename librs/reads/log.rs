// 读取log文件的脚本

use ndarray::Array2;
use regex::Regex;
use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
};

// 搜索标题行，找到每块数据的位置
pub fn search_title(path: &String) -> HashMap<String, u32> {
    let mut titles: HashMap<String, u32> = HashMap::new();

    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    for (l, line) in reader.lines().enumerate() {
        let line = line.unwrap();
        if line == "                          Input orientation:                          " {
            *titles.entry("geome".to_string()).or_insert(l as u32) = l as u32;
            continue;
        }
        if line == "                         Standard orientation:                         " {
            *titles.entry("geome".to_string()).or_insert(l as u32) = l as u32;
            continue;
        }
        if line == " AO basis set in the form of general basis input (Overlap normalization):" {
            *titles.entry("basis".to_string()).or_insert(l as u32) = l as u32;
            continue;
        }
        if line == " *** Overlap *** " {
            *titles.entry("smat".to_string()).or_insert(l as u32) = l as u32;
            continue;
        }
        if line == "     Molecular Orbital Coefficients:" {
            *titles.entry("cmat".to_string()).or_insert(l as u32) = l as u32;
            continue;
        }
        if line == "     Alpha Molecular Orbital Coefficients:" {
            *titles.entry("cmat_a".to_string()).or_insert(l as u32) = l as u32;
            continue;
        }
        if line == "     Beta Molecular Orbital Coefficients:" {
            *titles.entry("cmat_b".to_string()).or_insert(l as u32) = l as u32;
            continue;
        }
        let re = Regex::new(r" +NBasis = +\d+ +MinDer = \d+ +MaxDer = \d+").unwrap();
        if re.is_match(&line) {
            *titles.entry("nbas".to_string()).or_insert(l as u32) = l as u32;
            continue;
        }
    }
    println!("{:#?}", titles);
    return titles;
}

// 读取结构
pub fn read_geome(path: &String, title_number: u32) -> (Vec<u32>, Vec<[f64; 3]>) {
    let s1 = Regex::new(r" +\d+ +(\d+) +\d +(-?\d+.\d{6}) +(-?\d+.\d{6}) +(-?\d+.\d{6})").unwrap();
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let mut atms = Vec::new();
    let mut xyzs = Vec::new();
    for line in reader.lines().skip(title_number as usize + 5) {
        let line = line.unwrap();
        // println!("{}", line);
        if let Some(caps) = s1.captures(&line) {
            let atm: u32 = caps[1].parse().unwrap();
            let x: f64 = caps[2].parse().unwrap();
            let y: f64 = caps[3].parse().unwrap();
            let z: f64 = caps[4].parse().unwrap();
            atms.push(atm);
            xyzs.push([x, y, z]);
        } else {
            break;
        }
    }
    (atms, xyzs)
}

// 读取基函数数量
pub fn read_nbas(path: &String, title_number: u32) -> u32 {
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let mut nbas: u32 = 0;
    let re = Regex::new(r" +NBasis = +(\d+) +MinDer = \d+ +MaxDer = \d+").unwrap();
    for line in reader.lines().skip(title_number as usize) {
        let line = line.unwrap();
        // println!("{}", line);
        if let Some(caps) = re.captures(&line) {
            nbas = caps[1].parse::<u32>().unwrap();
        }
        break;
    }
    return nbas;
}

// 读取系数矩阵
pub fn read_cmat(
    path: &String,
    title_number: u32,
    nbas: u32,
) -> (
    Vec<u32>,
    Vec<u32>,
    Vec<String>,
    Vec<f64>,
    Vec<bool>,
    Array2<f64>,
) {
    let nblock = nbas / 5 + if nbas % 5 == 0 { 0 } else { 1 }; // 块的数量
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let mut cmat = Array2::<f64>::zeros((nbas as usize, nbas as usize));
    let mut ato_atms: Vec<u32> = Vec::with_capacity(nbas as usize);
    let mut ato_shls: Vec<u32> = Vec::with_capacity(nbas as usize);
    let mut ato_syms: Vec<String> = Vec::with_capacity(nbas as usize);
    let mut obt_engs: Vec<f64> = Vec::with_capacity(nbas as usize);
    let mut obt_occs: Vec<bool> = Vec::with_capacity(nbas as usize);
    // 先获取所有的行
    // println!("nblock: {nblock}");
    let mut lines = Vec::new();
    for (l, line) in reader.lines().skip(title_number as usize + 1).enumerate() {
        let line = line.unwrap();
        // println!("{l}: {line}");
        lines.push(line);
        if l > (nblock * (nbas + 3)) as usize {
            break;
        }
    }
    let mut iatm = 0;
    for i in 0..nblock {
        let loc = i * (nbas + 3) + 1;

        // 原子信息
        if i == 0 {
            for j in 0..nbas {
                let line = lines.get((loc + j + 2) as usize).unwrap();
                // println!("line:{}", line);
                let part = &line[4..15];
                // println!("part:{}", part);
                let atom_part = part[..4].trim();
                if !atom_part.is_empty() {
                    iatm = atom_part.parse::<u32>().unwrap();
                }
                ato_atms.push(iatm);

                let match_part = part[6..].trim();
                // println!("match_part:{}", match_part);
                let shl = match_part.chars().next().unwrap().to_string();
                let sym = match_part[1..].to_string();
                ato_shls.push(shl.parse::<u32>().unwrap());
                ato_syms.push(sym);
            }
        }

        // 轨道占据
        let occs_line = lines.get(loc as usize).unwrap();
        // println!("occs_line: {}", occs_line);
        let occs_part = &occs_line[21..];
        let occs: Vec<bool> = occs_part
            .chars()
            .collect::<Vec<_>>()
            .chunks(10)
            .map(|chunk| {
                let s: String = chunk.iter().collect();
                s.trim().ends_with('O')
            })
            .collect();
        obt_occs.extend_from_slice(&occs);

        // 轨道能量
        let engs_line = lines.get(loc as usize + 1).unwrap();
        // println!("engs_line: {}", engs_line);
        let engs_part = &engs_line[21..];
        let engs: Vec<f64> = engs_part
            .chars()
            .collect::<Vec<_>>()
            .chunks(10)
            .map(|chunk| {
                let s: String = chunk.iter().collect();
                s.trim().parse::<f64>().unwrap()
            })
            .collect();
        obt_engs.extend_from_slice(&engs);

        // 轨道系数
        for j in 0..nbas {
            let coef_line = lines.get((loc + j + 2) as usize).unwrap();
            let coef_part = &coef_line[21..];
            let coefs: Vec<f64> = coef_part
                .chars()
                .collect::<Vec<_>>()
                .chunks(10)
                .map(|chunk| {
                    let s: String = chunk.iter().collect();
                    // println!("{}", s);
                    s.trim().parse::<f64>().unwrap()
                })
                .collect();
            for (k, coef) in coefs.iter().enumerate() {
                cmat[(j as usize, i as usize * 5 + k)] = *coef;
            }
        }
    }
    (ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat)
}

pub fn read_basis(path: &String, title_number: u32) -> Vec<(u32, i32, u32, f64, f64)> {
    let s1 = Regex::new(r"^ +(\d+) +\d+$").unwrap();
    let s2 = Regex::new(r" ([SPD]+) +(\d+) \d.\d{2} +\d.\d{12}").unwrap();
    let s3 = Regex::new(r"^ +(( +-?\d.\d{10}D[+-]\d{2}){2,3})").unwrap();
    let mut ang_map: HashMap<char, u32> = HashMap::new();
    ang_map.insert('S', 0);
    ang_map.insert('P', 1);
    ang_map.insert('D', 2);
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let mut basis: Vec<(u32, i32, u32, f64, f64)> = Vec::new();
    let mut angs: Vec<u32> = Vec::new();
    let mut shl = 0;
    let mut atm: u32 = 0;
    for line in reader.lines().skip(title_number as usize + 1) {
        let line = line.unwrap();
        // println!("{}", line);
        if s1.is_match(&line) {
            // println!("{}", caps[1]);
            let caps = s1.captures(&line).unwrap();
            atm = caps[1].parse::<u32>().unwrap();
        } else if s2.is_match(&line) {
            let caps = s2.captures(&line).unwrap();
            let shl_sym: String = caps[1].to_string();
            angs = shl_sym
                .chars()
                .map(|s| *ang_map.get(&s).expect("未知角动量符号"))
                .collect();
            shl += 1;
        } else if s3.is_match(&line) {
            let caps = s3.captures(&line).unwrap();
            let nums: Vec<f64> = caps[1]
                .split_whitespace()
                .map(|each| {
                    let coef = each.replace("D", "e");
                    coef.parse::<f64>().unwrap()
                })
                .collect();
            if angs.len() == 1 {
                let alp = nums[0];
                let coe = nums[1];
                basis.push((atm, shl, angs[0], alp, coe));
            }
            if angs.len() == 2 {
                let alp = nums[0];
                let coe1 = nums[1];
                let coe2 = nums[2];
                basis.push((atm, shl, angs[0], alp, coe1));
                basis.push((atm, shl, angs[0], alp, coe2));
            }
        } else if line == " ****" {
            continue;
        } else {
            break;
        }
    }
    basis
}
