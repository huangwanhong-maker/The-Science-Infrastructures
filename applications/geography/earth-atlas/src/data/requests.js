// Data requested on demand from the local server, which has the laboratory's
// data library fetch a dataset with chosen parameters and build a product
// (data/acquisition/on_demand.py). The page asks, then follows the job until
// its product is ready; a page served without that server cannot ask.

export class Requests {
  constructor(base = '') {
    this.base = base;
    this.jobs = new Map();              // key -> { job, product, error, listeners }
  }

  // The state of a request by its key, if it has been made: { job, product, error }.
  get(key) {
    return this.jobs.get(key) || null;
  }

  // Asks for a dataset with parameters; calls onChange(entry) as the job advances.
  // Returns the entry; entry.product is set when the product has been read.
  ask(key, dataset, parameters, onChange) {
    const known = this.jobs.get(key);
    if (known && !known.error) return known;
    const entry = { job: null, product: null, error: null, asked: Date.now() };
    this.jobs.set(key, entry);
    const q = new URLSearchParams({ dataset, ...parameters });
    const follow = async (url) => {
      try {
        const response = await fetch(`${this.base}${url}`, { cache: 'no-store' });
        const job = await response.json().catch(() => ({ state: 'failed', message: `the server answered ${response.status}` }));
        entry.job = job;
        if (job.state === 'failed') throw new Error(job.message || 'the request failed');
        if (job.state === 'done') {
          const p = await fetch(`${this.base}data/${job.product}`, { cache: 'no-store' });
          if (!p.ok) throw new Error(`the product could not be read (${p.status})`);
          entry.product = await p.json();
          onChange(entry);
          return;
        }
        onChange(entry);
        setTimeout(() => follow(`request/${job.id}`), 3000);
      } catch (e) {
        entry.error = e.message === 'Failed to fetch' ? 'the local server could not be reached' : e.message;
        onChange(entry);
      }
    };
    follow(`request?${q}`);
    onChange(entry);
    return entry;
  }
}
