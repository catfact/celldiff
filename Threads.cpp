// create N number of threads
void init_threads() {
  // calculation counter lock
  // calculation-finished condition
  // copy-finished condition
}

// copy thread:
void * copy_thr() {
  // wait until calculations are complete (counter is at max)
  // copy data
  // wake up calculation threads
}


// calculate thread:
void * calculate_thr() {
  for (;;) {
    // wait on ready-to-compute condition
    // acquire, increment computation counter
  }
}
