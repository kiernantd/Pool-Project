use std::f64::consts::PI;

/// Basic stop (pool)
#[derive(Clone, Debug)]
struct Stop {
    id: usize,
    name: String,
    lat: f64, // degrees
    lon: f64, // degrees
    service_minutes: f64,
}

impl Stop {
    fn new(id: usize, name: &str, lat: f64, lon: f64, service_minutes: f64) -> Self {
        Self {
            id,
            name: name.to_string(),
            lat,
            lon,
            service_minutes,
        }
    }
}

/// Haversine distance in meters between two lat/lon points
fn haversine_meters(a_lat: f64, a_lon: f64, b_lat: f64, b_lon: f64) -> f64 {
    let r = 6_371_000.0; // earth radius in meters
    let to_rad = |deg: f64| deg * PI / 180.0;
    let (alat, alon, blat, blon) = (to_rad(a_lat), to_rad(a_lon), to_rad(b_lat), to_rad(b_lon));
    let dlat = blat - alat;
    let dlon = blon - alon;
    let sin_dlat = (dlat / 2.0).sin();
    let sin_dlon = (dlon / 2.0).sin();
    let a = sin_dlat * sin_dlat + alat.cos() * blat.cos() * sin_dlon * sin_dlon;
    let c = 2.0 * a.sqrt().asin();
    r * c
}

/// Convert meters to estimated drive minutes using an average speed (km/h)
fn meters_to_minutes(meters: f64, avg_kmph: f64) -> f64 {
    if avg_kmph <= 0.0 {
        return f64::INFINITY;
    }
    let km = meters / 1000.0;
    (km / avg_kmph) * 60.0
}

/// Route for a single cleaner
#[derive(Clone, Debug)]
struct Route {
    id: usize,
    stops: Vec<Stop>, // order matters
    depot: Option<Stop>, // optional start/end point (e.g., garage)
}

impl Route {
    fn new(id: usize, depot: Option<Stop>) -> Self {
        Self {
            id,
            stops: Vec::new(),
            depot,
        }
    }

    fn total_distance_meters(&self) -> f64 {
        let mut dist = 0.0;
        let mut prev_opt = self.depot.as_ref();
        for s in &self.stops {
            if let Some(prev) = prev_opt {
                dist += haversine_meters(prev.lat, prev.lon, s.lat, s.lon);
            }
            prev_opt = Some(s);
        }
        // return to depot if depot exists
        if let Some(depot) = &self.depot {
            if let Some(last) = self.stops.last() {
                dist += haversine_meters(last.lat, last.lon, depot.lat, depot.lon);
            }
        }
        dist
    }

    fn total_time_minutes(&self, avg_kmph: f64) -> f64 {
        let travel_meters = self.total_distance_meters();
        let travel_minutes = meters_to_minutes(travel_meters, avg_kmph);
        let service: f64 = self.stops.iter().map(|s| s.service_minutes).sum();
        travel_minutes + service
    }

    /// Nearest neighbor heuristic starting from depot (if present), else first stop
    fn build_nearest_neighbor(&mut self) {
        if self.stops.is_empty() {
            return;
        }
        let mut remaining = self.stops.clone();
        self.stops.clear();

        // determine starting point
        let mut current = if let Some(depot) = &self.depot {
            depot.clone()
        } else {
            remaining.remove(0)
        };

        // if depot was not originally in remaining (i.e., we used depot as start), ensure we don't include it
        loop {
            if remaining.is_empty() {
                break;
            }
            // find nearest stop to current
            let mut best_idx = 0usize;
            let mut best_dist = f64::INFINITY;
            for (i, s) in remaining.iter().enumerate() {
                let d = haversine_meters(current.lat, current.lon, s.lat, s.lon);
                if d < best_dist {
                    best_dist = d;
                    best_idx = i;
                }
            }
            let next = remaining.remove(best_idx);
            self.stops.push(next.clone());
            current = next;
        }
    }

    /// 2-opt improvement (simple implementation)
    fn two_opt(&mut self) {
        let n = self.stops.len();
        if n < 3 {
            return;
        }
        let mut improved = true;
        while improved {
            improved = false;
            for i in 0..n - 2 {
                for k in i + 2..n {
                    // avoid breaking adjacency with depot if needed (we consider route as path including depot as start and end)
                    if k == n - 1 && i == 0 && self.depot.is_none() {
                        continue;
                    }
                    let delta = self.two_opt_swap_delta(i, k);
                    if delta < -1e-6 {
                        self.do_two_opt_swap(i, k);
                        improved = true;
                    }
                }
            }
        }
    }

    /// compute change in distance if we reverse segment (i+1..=k)
    fn two_opt_swap_delta(&self, i: usize, k: usize) -> f64 {
        // nodes: A - B ... C - D
        // edges removed: AB and CD
        // edges added: AC and BD
        let n = self.stops.len();
        let get_point = |idx_opt: Option<isize>| -> (f64, f64) {
            match idx_opt {
                Some(idx) if idx >= 0 && (idx as usize) < n => {
                    let s = &self.stops[idx as usize];
                    (s.lat, s.lon)
                }
                _ => {
                    // depot or out-of-range -> use depot (if present) otherwise repeat endpoint
                    if let Some(depot) = &self.depot {
                        (depot.lat, depot.lon)
                    } else if n == 0 {
                        (0.0, 0.0)
                    } else {
                        let s = &self.stops[n - 1];
                        (s.lat, s.lon)
                    }
                }
            }
        };

        let a_idx = if i == 0 {
            None
        } else {
            Some((i - 1) as isize)
        };
        let b_idx = Some(i as isize);
        let c_idx = Some(k as isize);
        let d_idx = if k + 1 >= n { None } else { Some((k + 1) as isize) };

        let (a_lat, a_lon) = get_point(a_idx);
        let (b_lat, b_lon) = get_point(b_idx);
        let (c_lat, c_lon) = get_point(c_idx);
        let (d_lat, d_lon) = get_point(d_idx);

        let removed = haversine_meters(a_lat, a_lon, b_lat, b_lon) + haversine_meters(c_lat, c_lon, d_lat, d_lon);
        let added = haversine_meters(a_lat, a_lon, c_lat, c_lon) + haversine_meters(b_lat, b_lon, d_lat, d_lon);
        added - removed
    }

    fn do_two_opt_swap(&mut self, i: usize, k: usize) {
        // reverse the segment i..=k in place
        self.stops[i..=k].reverse();
    }

    /// Convenience: optimize by building NN then 2-opt
    fn optimize(&mut self) {
        if self.stops.is_empty() {
            return;
        }
        // copy current to temp and run NN on that set
        let all_stops = self.stops.clone();
        self.stops = all_stops;
        self.build_nearest_neighbor();
        self.two_opt();
    }

    fn print_summary(&self, avg_kmph: f64) {
        println!("Route {} summary:", self.id);
        println!("  Stops ({}):", self.stops.len());
        for s in &self.stops {
            println!("    {}: {} ({:.6}, {:.6})", s.id, s.name, s.lat, s.lon);
        }
        let meters = self.total_distance_meters();
        let minutes = self.total_time_minutes(avg_kmph);
        println!("  Total distance: {:.1} m", meters);
        println!("  Estimated time (incl service): {:.1} min (avg {:.1} km/h)", minutes, avg_kmph);
    }
}

/// Manager for multiple routes (cleaners)
struct RouteManager {
    routes: Vec<Route>,
}

impl RouteManager {
    fn new() -> Self {
        Self { routes: Vec::new() }
    }

    fn add_route(&mut self, route: Route) {
        self.routes.push(route);
    }

    /// Add stop to a specific route (by route id)
    fn add_stop_to_route(&mut self, route_id: usize, stop: Stop) {
        if let Some(r) = self.routes.iter_mut().find(|r| r.id == route_id) {
            r.stops.push(stop);
        } else {
            eprintln!("Route {} not found", route_id);
        }
    }

    /// Remove stop by stop id from any route (returns the stop if found)
    fn remove_stop_by_id(&mut self, stop_id: usize) -> Option<Stop> {
        for r in self.routes.iter_mut() {
            if let Some(pos) = r.stops.iter().position(|s| s.id == stop_id) {
                return Some(r.stops.remove(pos));
            }
        }
        None
    }

    /// Reassign stop to another route by ids
    fn reassign_stop(&mut self, stop_id: usize, target_route_id: usize) -> bool {
        if let Some(stop) = self.remove_stop_by_id(stop_id) {
            self.add_stop_to_route(target_route_id, stop);
            true
        } else {
            false
        }
    }

    /// Optimize all routes
    fn optimize_all(&mut self) {
        for r in self.routes.iter_mut() {
            r.optimize();
        }
    }

    fn print_all_summaries(&self, avg_kmph: f64) {
        for r in &self.routes {
            r.print_summary(avg_kmph);
            println!();
        }
    }
}

/// Example usage: create some synthetic stops and demonstrate add/remove/reassign/optimize.
fn main() {
    // average speed assumption for drive time estimates
    let avg_kmph = 35.0;

    // create a (synthetic) depot (e.g., company yard)
    let depot = Stop::new(0, "Depot", 40.4406, -79.9959, 0.0); // Pittsburgh-ish coords

    // sample stops (lat, lon roughly in Pittsburgh area). In real use, read from CSV or API.
    let sample_stops = vec![
        Stop::new(1, "Pool A", 40.4475, -79.9646, 10.0),
        Stop::new(2, "Pool B", 40.4300, -80.0005, 12.0),
        Stop::new(3, "Pool C", 40.4305, -79.9800, 8.0),
        Stop::new(4, "Pool D", 40.4520, -79.9730, 15.0),
        Stop::new(5, "Pool E", 40.4200, -79.9800, 10.0),
        Stop::new(6, "Pool F", 40.4380, -80.0100, 10.0),
        Stop::new(7, "Pool G", 40.4450, -80.0050, 10.0),
    ];

    // create two routes/cleaners
    let mut manager = RouteManager::new();
    let mut r1 = Route::new(1, Some(depot.clone()));
    let mut r2 = Route::new(2, Some(depot.clone()));

    // naive split: first N/2 go to route1, rest to route2
    for (i, s) in sample_stops.into_iter().enumerate() {
        if i % 2 == 0 {
            r1.stops.push(s);
        } else {
            r2.stops.push(s);
        }
    }

    manager.add_route(r1);
    manager.add_route(r2);

    println!("Before optimization:");
    manager.print_all_summaries(avg_kmph);

    // Optimize all routes
    manager.optimize_all();

    println!("After optimization:");
    manager.print_all_summaries(avg_kmph);

    // Demonstrate add/remove/reassign:
    println!("Demonstrating add/remove/reassign:");

    // Add a new stop to route 1
    let new_stop = Stop::new(99, "Pool Z", 40.4350, -79.9750, 10.0);
    manager.add_stop_to_route(1, new_stop);
    manager.routes.iter_mut().find(|r| r.id == 1).unwrap().optimize();

    println!("[After adding Pool Z to Route 1 and optimizing]");
    manager.print_all_summaries(avg_kmph);

    // Remove stop id 3
    if let Some(removed) = manager.remove_stop_by_id(3) {
        println!("Removed stop {} ({})", removed.id, removed.name);
    }

    // Reassign stop id 2 to route 1
    if manager.reassign_stop(2, 1) {
        println!("Reassigned stop 2 to Route 1");
    }

    // Re-optimize
    manager.optimize_all();
    println!("[Final optimized routes]");
    manager.print_all_summaries(avg_kmph);
}
