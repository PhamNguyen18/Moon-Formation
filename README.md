# Moon-Formation 
THIS REPOSITORY IS NO LONGER MAINTAINED

A work in progress lunar accretion simulation. This code takes inspiration from the work of Salmon and Canup 2012 who used a hybrid hydro/ N-body code to model the accretion of the Moon. A 1D fluid is used to model the disk within the Roche limit and an N-body code outside of it. For this project I am using the open source N-body code [rebound](https://github.com/hannorein/rebound). The hydrocode will be custom built at a later time. 

## Recently added: Experimental merging routine that includes tidal forces

Still trouble shooting this new routine as it currently does not work. The routine can be called by adding the following to your problem file:

r->collision_resolve = reb_collision_resolve_tidal;

TODO:
- Replicate results of Ida 1997
- Implement hydrocode
- Add drag forces to planetesimals 
- Come up with a clever repo name 
