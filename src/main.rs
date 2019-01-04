#[macro_use] extern crate lazy_static;

use std::ops::{Add, Sub, Mul, Rem};
use std::io::Write;
use std::fs;

fn random() -> f64 {
    rand::random::<f64>()
}

#[derive(Copy, Clone, Debug)]
struct Vec {
    x: f64,
    y: f64,
    z: f64,
}

impl Vec {
    fn new(x: f64, y: f64, z: f64) -> Vec {
        Vec {x, y, z}
    }
    fn zero() -> Vec {
        Vec::new(0.0, 0.0, 0.0)
    }
    fn mult(&self, b: &Vec) -> Vec {
        Vec::new(self.x * b.x, self.y * b.y, self.z * b.z)
    }
    fn norm(mut self) -> Vec {
        let l = 1.0 / (self.x * self.x + self.y * self.y + self.z * self.z).sqrt();
        self.x = self.x * l;
        self.y = self.y * l;
        self.z = self.z * l;
        self
    }
    fn dot(&self, b: &Vec) -> f64 {
        return self.x * b.x + self.y * b.y + self.z * b.z;
    }
}

impl Add for Vec {
    type Output = Vec;
    fn add(self, rhs: Self) -> Self {
        Vec::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl Sub for Vec {
    type Output = Vec;
    fn sub(self, rhs: Self) -> Self {
        Vec::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl Mul<f64> for Vec {
    type Output = Vec;
    fn mul(self, rhs: f64) -> Self {
        Vec::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl Rem for Vec {
    type Output = Vec;
    fn rem(self, rhs: Self) -> Self {
        Vec::new(
            self.y * rhs.z - self.z * rhs.y,
            self.z * rhs.x - self.x * rhs.z,
            self.x * rhs.y - self.y * rhs.x
        )
    }
}

#[derive(Debug)]
struct Ray {
    o: Vec,
    d: Vec,
}

impl Ray {
    fn new(o: Vec, d: Vec) -> Ray {
        Ray { o, d }
    }
}

enum Refl {
    Diff,
    Spec,
    Refr,
}

struct Sphere {
    rad: f64,
    p: Vec,
    e: Vec,
    c: Vec,
    refl: Refl,
}

impl Sphere {
    fn intersect(&self, r: &Ray) -> f64 {
        let op = self.p - r.o;
        let eps = 1e-4;
        let b = op.dot(&r.d);
        let mut det = b * b - op.dot(&op) + self.rad * self.rad;
        if det < 0.0 {
            return 0.0;
        }
        det = det.sqrt();
        let t = b - det;
        if t > eps {
            return t;
        }
        let t = b + det;
        if t > eps {
            return t;
        } else {
            0.0
        }
    }
}

lazy_static! {
    static ref spheres: [Sphere; 9] = [
        Sphere { rad: 1e5, p: Vec::new(1e5 + 1.0, 40.8, 81.6), e: Vec::zero(), c: Vec::new(0.75, 0.25, 0.25), refl: Refl::Diff },
        Sphere { rad: 1e5, p: Vec::new(-1e5 + 99.0, 40.8, 81.6), e: Vec::zero(), c: Vec::new(0.25, 0.25, 0.75), refl: Refl::Diff },
        Sphere { rad: 1e5, p: Vec::new(50.0, 40.8, 1e5), e: Vec::zero(), c: Vec::new(0.75, 0.75, 0.75), refl: Refl::Diff },
        Sphere { rad: 1e5, p: Vec::new(50.0, 40.8, -1e5 + 170.0), e: Vec::zero(), c: Vec::zero(), refl: Refl::Diff },
        Sphere { rad: 1e5, p: Vec::new(50.0, 1e5, 81.6), e: Vec::zero(), c: Vec::new(0.75, 0.75, 0.75), refl: Refl::Diff },
        Sphere { rad: 1e5, p: Vec::new(50.0, -1e5 + 81.6, 81.6), e: Vec::zero(), c: Vec::new(0.75, 0.75, 0.75), refl: Refl::Diff },
        Sphere { rad: 16.5, p: Vec::new(27.0, 16.5, 47.0), e: Vec::zero(), c: Vec::new(1.0, 1.0, 1.0) * 0.999, refl: Refl::Spec },
        Sphere { rad: 16.5, p: Vec::new(73.0, 16.5, 78.0), e: Vec::zero(), c: Vec::new(1.0, 1.0, 1.0) * 0.999, refl: Refl::Refr },
        Sphere { rad: 600.0, p: Vec::new(50.0, 681.6 - 0.27, 81.6), e: Vec::new(12.0, 12.0, 12.0), c: Vec::zero(), refl: Refl::Diff },
    ];
}

fn clamp(x: f64) -> f64 {
    if x < 0.0 {
        0.0
    } else if x > 1.0 {
        1.0
    } else {
        x
    }
}

fn to_int(x: f64) -> u8 {
    (f64::powf(clamp(x), 1.0 / 2.2) * 255.0 + 0.5) as u8
}

fn intersect(r: &Ray, t: &mut f64, id: &mut usize) -> bool {
    let n = spheres.len();
    let inf: f64 = 1e20;
    *t = inf;
    for i in (0..n).rev() {
        let d = spheres[i].intersect(r);
        if d != 0.0 && d < *t {
            *t = d;
            *id = i;
        }
    }
    return *t < inf;
}

fn radiance(r: &Ray, depth: u8) -> Vec {
    let mut t: f64 = 0.0;
    let mut id = 0;
    if !intersect(r, &mut t, &mut id) {
        return Vec::zero();
    }
    let obj = &spheres[id];
    let x = r.o + r.d * t;
    let n = (x - obj.p).norm();
    let nl = if n.dot(&r.d) < 0.0 { n } else { n * -1.0 };
    let mut f = obj.c;
    let p = if f.x > f.y && f.x > f.z { f.x } else if f.y > f.z { f.y } else { f.z };
    let depth = depth + 1;
    if depth > 5 {
        if depth < 127 && random() < p {
            f = f * (1.0 / p);
        } else {
            return obj.e;
        }
    }

    return match obj.refl {
        Refl::Diff => {
            let r1 = 2.0 * std::f64::consts::PI * random();
            let r2 = random();
            let r2s = r2.sqrt();
            let w = nl;
            let u = ((if w.x.abs() > 0.1 { Vec::new(0.0, 1.0, 0.0) } else { Vec::new(1.0, 0.0, 0.0) }) % w).norm();
            let v = w % u;
            let d = (u * f64::cos(r1) * r2s + v * f64::sin(r1) * r2s + w * (1.0 - r2).sqrt()).norm();
            obj.e + f.mult(&radiance(&Ray::new(x, d), depth))
        },
        Refl::Spec => {
            obj.e + f.mult(&radiance(&Ray::new(x, r.d - n * 2.0 * n.dot(&r.d)), depth))
        },
        _ => { // Refl.Refr
            let refl_ray = Ray::new(x, r.d - n * 2.0 * n.dot(&r.d));
            let into = n.dot(&nl) > 0.0;
            let nc = 1.0;
            let nt = 1.5;
            let nnt = if into { nc / nt } else { nt / nc };
            let ddn = r.d.dot(&nl);
            let cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);
            if cos2t < 0.0 {
                obj.e + f.mult(&radiance(&refl_ray, depth))
            } else {
                let tdir = r.d * nnt - n * ((if into { 1.0 } else { -1.0 }) * (ddn * nnt + cos2t.sqrt()));
                tdir.norm();
                let a = nt - nc;
                let b = nt + nc;
                let r0 = a * a / (b * b);
                let c = 1.0 - (if into { -ddn } else {tdir.dot(&n)});
                let re = r0 + (1.0 - r0) * c * c * c * c * c;
                let tr = 1.0 - re;
                let p = 0.25 + 0.5 * re;
                let rp = re / p;
                let tp = tr / (1.0 - p);
                obj.e + f.mult(&(
                    if depth > 2 {
                        if random() < p {
                            radiance(&refl_ray, depth) * rp
                        } else {
                            radiance(&Ray::new(x, tdir), depth) * tp
                        }
                    } else {
                        radiance(&refl_ray, depth) * re + radiance(&Ray::new(x, tdir), depth) * tr
                    }
                ))
            }
        }
    }

}

fn main() {
    let w = 1024;
    let h = 768;
    let samps = if std::env::args().len() == 2 { std::env::args().skip(1).next().unwrap().parse().unwrap() } else { 1 };
    let cam = Ray::new(Vec::new(50.0, 52.0, 295.6), Vec::new(0.0, -0.042612, -1.0).norm());

    let cx = Vec::new((w as f64) * 0.5135 / (h as f64), 0.0, 0.0);
    let cy = (cx % cam.d).norm() * 0.5135;
    let mut c = vec![Vec::zero(); w * h];

    for y in 0..h {
        writeln!(std::io::stderr(), "Rendering ({} spp) {:5.2}%", samps * 4, 100.0 * (y as f64) / ((h as f64) - 1.0)).unwrap();
        for x in 0..w {
            let i = (h - y - 1) * w + x;
            let mut r = Vec::zero();
            for sy in 0..2 {
                for sx in 0..2 {
                    for _s in 0..samps {
                        let r1 = 2.0 * random();
                        let dx = if r1 < 1.0 { r1.sqrt() - 1.0 } else { 1.0 - (2.0 - r1).sqrt() };
                        let r2 = 2.0 * random();
                        let dy = if r2 < 1.0 { r2.sqrt() - 1.0 } else { 1.0 - (2.0 - r2).sqrt() }; 
                        let d = cx * ((((sx as f64) + 0.5 + dx) / 2.0 + (x as f64)) / (w as f64) - 0.5)
                            + cy * ((((sy as f64) + 0.5 + dy) / 2.0 + (y as f64)) / (h as f64) - 0.5) + cam.d;
                        r = r + radiance(&(Ray::new(cam.o + d * 140.0, d.norm())), 0) * (1.0 / (samps as f64));
                    }
                    c[i] = c[i] + Vec::new(clamp(r.x), clamp(r.y), clamp(r.z)) * 0.25;
                }
            }
        }
    }


    let mut f = fs::File::create("image.ppm").unwrap();
    writeln!(f, "P3\n{} {}\n{}", w, h, 255).unwrap();
    for i in 0..w*h {
        writeln!(f, "{} {} {}", to_int(c[i].x), to_int(c[i].y), to_int(c[i].z)).unwrap();
    }
}
