import { useState, useEffect, createContext, useContext } from 'react';

/**
 * Context for database access throughout the app.
 */
export const DatabaseContext = createContext(null);

/**
 * Hook to access the database context.
 */
export function useDatabase() {
  const context = useContext(DatabaseContext);
  if (!context) {
    throw new Error('useDatabase must be used within DatabaseProvider');
  }
  return context;
}

/**
 * Database URL configuration.
 * In production, fetches from the hosted site.
 * Falls back to local path for development or if production URL fails.
 */
const PRODUCTION_DATABASE_URL = 'https://waudbylab.org/nmr-pH/database/current/database.json';
const LOCAL_DATABASE_URL = './database/current/database.json';

// Use production URL by default, with local fallback
const DEFAULT_DATABASE_URL = PRODUCTION_DATABASE_URL;

/**
 * DatabaseLoader component.
 * Fetches the buffer database and provides it via context.
 */
export function DatabaseLoader({ children, databaseUrl = DEFAULT_DATABASE_URL }) {
  const [database, setDatabase] = useState(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);

  useEffect(() => {
    let cancelled = false;

    async function loadDatabase() {
      try {
        setLoading(true);
        setError(null);

        let data;
        let fetchError = null;

        // Try primary URL first
        try {
          const response = await fetch(databaseUrl);
          if (!response.ok) {
            throw new Error(`Failed to load database: ${response.status} ${response.statusText}`);
          }
          data = await response.json();
        } catch (err) {
          fetchError = err;
          // If primary URL fails and we're using production URL, try local fallback
          if (databaseUrl === PRODUCTION_DATABASE_URL) {
            console.warn('Production database fetch failed, trying local fallback:', err);
            try {
              const response = await fetch(LOCAL_DATABASE_URL);
              if (!response.ok) {
                throw new Error(`Failed to load local database: ${response.status} ${response.statusText}`);
              }
              data = await response.json();
              fetchError = null; // Clear error since local fetch succeeded
            } catch (localErr) {
              console.warn('Local database fetch also failed:', localErr);
              throw fetchError; // Throw original error
            }
          } else {
            throw err;
          }
        }

        if (!cancelled) {
          // Build maps for quick lookup
          const samplesMap = new Map(data.samples.map(s => [s.sample_id, s]));
          const buffersMap = new Map(data.buffers.map(b => [b.buffer_id, b]));

          setDatabase({
            ...data,
            samplesMap,
            buffersMap
          });
          setLoading(false);

          // Cache in localStorage for offline use
          try {
            localStorage.setItem('nmr-ph-database', JSON.stringify(data));
            localStorage.setItem('nmr-ph-database-timestamp', Date.now().toString());
          } catch (e) {
            console.warn('Failed to cache database:', e);
          }
        }
      } catch (err) {
        console.error('Database load error:', err);

        // Try to load from cache
        try {
          const cached = localStorage.getItem('nmr-ph-database');
          if (cached) {
            const data = JSON.parse(cached);
            const samplesMap = new Map(data.samples.map(s => [s.sample_id, s]));
            const buffersMap = new Map(data.buffers.map(b => [b.buffer_id, b]));

            if (!cancelled) {
              setDatabase({
                ...data,
                samplesMap,
                buffersMap
              });
              setLoading(false);
              setError('Using cached database (could not fetch latest)');
            }
            return;
          }
        } catch (cacheErr) {
          console.warn('Cache load failed:', cacheErr);
        }

        if (!cancelled) {
          setError(err.message);
          setLoading(false);
        }
      }
    }

    loadDatabase();

    return () => {
      cancelled = true;
    };
  }, [databaseUrl]);

  // Get available solvents from database
  const solvents = database
    ? [...new Set(database.samples.map(s => s.solvent))].filter(Boolean)
    : [];

  // Get buffers filtered by solvent
  const getBuffersForSolvent = (solvent) => {
    if (!database) return [];
    const sampleIds = database.samples
      .filter(s => s.solvent === solvent)
      .map(s => s.sample_id);
    return database.buffers.filter(b => sampleIds.includes(b.sample_id));
  };

  // Get available nuclei for selected buffers
  const getNucleiForBuffers = (buffers) => {
    const nuclei = new Set();
    for (const buffer of buffers) {
      for (const nucleus of Object.keys(buffer.chemical_shifts)) {
        nuclei.add(nucleus);
      }
    }
    return [...nuclei].sort();
  };

  const contextValue = {
    database,
    loading,
    error,
    solvents,
    getBuffersForSolvent,
    getNucleiForBuffers,
    getSample: (sampleId) => database?.samplesMap.get(sampleId),
    getBuffer: (bufferId) => database?.buffersMap.get(bufferId)
  };

  if (loading) {
    return (
      <div className="database-loader loading">
        <div className="spinner"></div>
        <p>Loading buffer database...</p>
      </div>
    );
  }

  if (error && !database) {
    return (
      <div className="database-loader error">
        <h2>Database Load Error</h2>
        <p>{error}</p>
        <p>
          Please check your internet connection or{' '}
          <a href="https://github.com/waudbylab/nmr-pH/issues" target="_blank" rel="noopener noreferrer">
            report an issue
          </a>.
        </p>
      </div>
    );
  }

  return (
    <DatabaseContext.Provider value={contextValue}>
      {error && (
        <div className="database-warning">
          <span className="warning-icon">&#9888;</span> {error}
        </div>
      )}
      {children}
    </DatabaseContext.Provider>
  );
}

export default DatabaseLoader;
