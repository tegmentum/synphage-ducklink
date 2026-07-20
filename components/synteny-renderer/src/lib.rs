//! bio-synteny-renderer -- SVG synteny plots as a DuckLink-loadable
//! component AND a reusable rendering capability.
//!
//! Exports two WIT surfaces from one .wasm:
//! - `tegmentum:bio/synteny-renderer` -- the biological capability,
//!   callable by any component-model host that wants tracks/features/links
//!   in and SVG bytes out.
//! - `duckdb:extension/{guest, callback-dispatch, table-stream-dispatch}`
//!   -- DuckLink's dispatch surface. One table function,
//!   `render_synteny_svg`, is registered during `load()` and delegates to
//!   the same underlying `render::render_svg` implementation.
//!
//! Both the biology Guest impl and the DuckLink bridge share one
//! `Component` type; everything interesting lives in `render.rs`. The two
//! surfaces are cfg-gated to `target_arch = "wasm32"` so `cargo test` on
//! the host target keeps exercising the pure `render` module without
//! dragging in the wasm-only bindings.

#[cfg(target_arch = "wasm32")]
#[allow(warnings)]
mod bindings;

pub mod render;

#[cfg(target_arch = "wasm32")]
mod ducklink;

#[cfg(target_arch = "wasm32")]
pub(crate) struct Component;

#[cfg(target_arch = "wasm32")]
mod biology {
    use crate::bindings::exports::tegmentum::bio::synteny_renderer::{
        Feature as WitFeature, Guest, Link as WitLink, RenderError as WitRenderError,
        Track as WitTrack,
    };
    use crate::render::{self, Feature, Link, RenderError, Track};
    use crate::Component;

    impl Guest for Component {
        fn render_svg(
            tracks: Vec<WitTrack>,
            features: Vec<WitFeature>,
            links: Vec<WitLink>,
        ) -> Result<Vec<u8>, WitRenderError> {
            let tracks: Vec<Track> = tracks.into_iter().map(track_from_wit).collect();
            let features: Vec<Feature> = features.into_iter().map(feature_from_wit).collect();
            let links: Vec<Link> = links.into_iter().map(link_from_wit).collect();

            render::render_svg(&tracks, &features, &links).map_err(error_to_wit)
        }
    }

    fn track_from_wit(t: WitTrack) -> Track {
        Track {
            track_id: t.track_id,
            label: t.label,
            length: t.length,
        }
    }

    fn feature_from_wit(f: WitFeature) -> Feature {
        Feature {
            track_id: f.track_id,
            feature_id: f.feature_id,
            start_position: f.start_position,
            end_position: f.end_position,
            strand: f.strand,
            colour: f.colour,
            label: f.label,
        }
    }

    fn link_from_wit(l: WitLink) -> Link {
        Link {
            query_track: l.query_track,
            query_feature: l.query_feature,
            subject_track: l.subject_track,
            subject_feature: l.subject_feature,
            identity: l.identity,
            colour: l.colour,
        }
    }

    fn error_to_wit(e: RenderError) -> WitRenderError {
        match e {
            RenderError::InvalidModel(msg) => WitRenderError::InvalidModel(msg),
            RenderError::EmptyInput => WitRenderError::EmptyInput,
        }
    }
}

#[cfg(target_arch = "wasm32")]
bindings::export!(Component with_types_in bindings);
