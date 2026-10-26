#ifndef PhysicsTools_SimpleVertexTable_TwoTrackVertexFitter_h
#define PhysicsTools_SimpleVertexTable_TwoTrackVertexFitter_h

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <limits>
#include <vector>

#include "TMath.h"

// Reusable, stream-local two-track vertex fit and derived-geometry helper.
// Position and covariance use the standard CMSSW units: cm and cm^2.
class TwoTrackVertexFitter {
public:
  static constexpr float kInvalid = -9999.f;

  struct Result {
    int status = 0;  // 1 success, 2 invalid track state, 3 invalid fitter result, 4 exception/invalid numerics
    int fitValid = 0;
    int pvValid = 0;
    int lxyValid = 0;
    int signedLxyValid = 0;
    int l3dValid = 0;
    int pointing2DValid = 0;
    int pointing3DValid = 0;
    int perpendicularValid = 0;
    float chi2 = kInvalid, ndof = kInvalid, chi2Ndof = kInvalid, probability = kInvalid;
    float x = kInvalid, y = kInvalid, z = kInvalid;
    std::array<float, 6> covariance{{kInvalid, kInvalid, kInvalid, kInvalid, kInvalid, kInvalid}};
    float sigmaX = kInvalid, sigmaY = kInvalid, sigmaZ = kInvalid;
    float lxy = kInvalid, lxyError = kInvalid, lxySignificance = kInvalid;
    float signedLxy = kInvalid, signedLxySignificance = kInvalid;
    float l3d = kInvalid, l3dError = kInvalid, l3dSignificance = kInvalid;
    float cosPointing2D = kInvalid, cosPointing3D = kInvalid;
    float perpendicular = kInvalid, perpendicularSignificance = kInvalid;
  };

  TwoTrackVertexFitter() : fitter_(false) {}

  Result fit(const reco::TransientTrack& first,
             const reco::TransientTrack& second,
             const reco::Track& firstTrack,
             const reco::Track& secondTrack,
             const reco::Vertex* primaryVertex,
             std::vector<reco::TransientTrack>& scratchTracks) const {
    Result out;
    if (!first.isValid() || !second.isValid()) {
      out.status = 2;
      return out;
    }

    try {
      // The caller reserves two slots once per event and reuses this container.
      scratchTracks.clear();
      scratchTracks.push_back(first);
      scratchTracks.push_back(second);
      const reco::Vertex fitted = fitter_.vertex(scratchTracks);
      if (!fitted.isValid()) {
        out.status = 3;
        return out;
      }

      const double chi2 = fitted.chi2();
      const double ndof = fitted.ndof();
      if (!std::isfinite(chi2) || !std::isfinite(ndof) || ndof > std::numeric_limits<int>::max()) {
        out.status = 4;
        return out;
      }
      const int ndofInt = static_cast<int>(std::lround(ndof));
      const auto& cov = fitted.covariance();
      const double values[] = {fitted.x(), fitted.y(), fitted.z(), chi2, ndof,
                               cov(0, 0), cov(0, 1), cov(0, 2), cov(1, 1), cov(1, 2), cov(2, 2)};
      for (double value : values) {
        if (!std::isfinite(value)) {
          out.status = 4;
          return out;
        }
      }
      if (ndofInt <= 0 || std::abs(ndof - ndofInt) > 1.e-6 || chi2 < 0. ||
          cov(0, 0) < 0. || cov(1, 1) < 0. || cov(2, 2) < 0.) {
        out.status = 4;
        return out;
      }

      out.status = 1;
      out.fitValid = 1;
      out.x = fitted.x(); out.y = fitted.y(); out.z = fitted.z();
      out.chi2 = chi2; out.ndof = ndof;
      out.chi2Ndof = chi2 / ndof;
      out.probability = TMath::Prob(chi2, ndofInt);
      if (!std::isfinite(out.chi2Ndof) || !std::isfinite(out.probability) || out.probability < 0.f || out.probability > 1.f) {
        out = Result{};
        out.status = 4;
        return out;
      }
      out.covariance = {{static_cast<float>(cov(0, 0)), static_cast<float>(cov(0, 1)),
                         static_cast<float>(cov(0, 2)), static_cast<float>(cov(1, 1)),
                         static_cast<float>(cov(1, 2)), static_cast<float>(cov(2, 2))}};
      out.sigmaX = std::sqrt(cov(0, 0)); out.sigmaY = std::sqrt(cov(1, 1)); out.sigmaZ = std::sqrt(cov(2, 2));
      const float rawFloats[] = {out.x, out.y, out.z, out.chi2, out.ndof, out.chi2Ndof, out.probability,
                                 out.sigmaX, out.sigmaY, out.sigmaZ, out.covariance[0], out.covariance[1],
                                 out.covariance[2], out.covariance[3], out.covariance[4], out.covariance[5]};
      for (float value : rawFloats) {
        if (!std::isfinite(value)) {
          out = Result{};
          out.status = 4;
          return out;
        }
      }

      if (!primaryVertex || primaryVertex->isFake()) return out;
      const auto& pvc = primaryVertex->covariance();
      const double dx = out.x - primaryVertex->x();
      const double dy = out.y - primaryVertex->y();
      const double dz = out.z - primaryVertex->z();
      const double lxy = std::hypot(dx, dy);
      const double l3d = std::sqrt(dx * dx + dy * dy + dz * dz);
      const double px = firstTrack.px() + secondTrack.px();
      const double py = firstTrack.py() + secondTrack.py();
      const double pz = firstTrack.pz() + secondTrack.pz();
      const double pt = std::hypot(px, py);
      const double p = std::sqrt(px * px + py * py + pz * pz);
      const double cxx = cov(0, 0) + pvc(0, 0);
      const double cxy = cov(0, 1) + pvc(0, 1);
      const double cxz = cov(0, 2) + pvc(0, 2);
      const double cyy = cov(1, 1) + pvc(1, 1);
      const double cyz = cov(1, 2) + pvc(1, 2);
      const double czz = cov(2, 2) + pvc(2, 2);
      const double lxyVar = lxy > 0. ? (dx * dx * cxx + 2. * dx * dy * cxy + dy * dy * cyy) / (lxy * lxy) : -1.;
      const double l3dVar = l3d > 0. ? (dx * dx * cxx + 2. * dx * dy * cxy + 2. * dx * dz * cxz +
                                           dy * dy * cyy + 2. * dy * dz * cyz + dz * dz * czz) / (l3d * l3d) : -1.;
      if (!std::isfinite(dx + dy + dz + lxy + l3d + cxx + cxy + cxz + cyy + cyz + czz)) return out;
      out.pvValid = 1;

      if (lxy > 0. && lxyVar > 0. && std::isfinite(lxyVar)) {
        const double error = std::sqrt(lxyVar), significance = lxy / error;
        if (std::isfinite(error) && std::isfinite(significance) && std::isfinite(static_cast<float>(lxy)) &&
            std::isfinite(static_cast<float>(error)) && std::isfinite(static_cast<float>(significance))) {
          out.lxy = lxy; out.lxyError = error; out.lxySignificance = significance; out.lxyValid = 1;
        }
      }
      if (l3d > 0. && l3dVar > 0. && std::isfinite(l3dVar)) {
        const double error = std::sqrt(l3dVar), significance = l3d / error;
        if (std::isfinite(error) && std::isfinite(significance) && std::isfinite(static_cast<float>(l3d)) &&
            std::isfinite(static_cast<float>(error)) && std::isfinite(static_cast<float>(significance))) {
          out.l3d = l3d; out.l3dError = error; out.l3dSignificance = significance; out.l3dValid = 1;
        }
      }
      if (pt > 0. && std::isfinite(pt)) {
        const double ux = px / pt, uy = py / pt;
        const double signedLxy = dx * ux + dy * uy;
        const double alongVar = ux * ux * cxx + 2. * ux * uy * cxy + uy * uy * cyy;
        if (alongVar > 0. && std::isfinite(alongVar)) {
          const double significance = signedLxy / std::sqrt(alongVar);
          if (std::isfinite(signedLxy) && std::isfinite(significance) && std::isfinite(static_cast<float>(signedLxy)) &&
              std::isfinite(static_cast<float>(significance))) {
            out.signedLxy = signedLxy;
            out.signedLxySignificance = significance;
            out.signedLxyValid = 1;
          }
        }
        if (lxy > 0.) {
          const double cosine = std::clamp((dx * px + dy * py) / (lxy * pt), -1., 1.);
          if (std::isfinite(cosine)) { out.cosPointing2D = cosine; out.pointing2DValid = 1; }
        }
        const double nx = -uy, ny = ux;
        const double perpVar = nx * nx * cxx + 2. * nx * ny * cxy + ny * ny * cyy;
        if (perpVar > 0. && std::isfinite(perpVar)) {
          const double value = dx * nx + dy * ny;
          const double significance = value / std::sqrt(perpVar);
          if (std::isfinite(value) && std::isfinite(significance) && std::isfinite(static_cast<float>(value)) &&
              std::isfinite(static_cast<float>(significance))) {
            out.perpendicular = value;
            out.perpendicularSignificance = significance;
            out.perpendicularValid = 1;
          }
        }
      }
      if (p > 0. && l3d > 0. && std::isfinite(p)) {
        const double cosine = std::clamp((dx * px + dy * py + dz * pz) / (l3d * p), -1., 1.);
        if (std::isfinite(cosine)) { out.cosPointing3D = cosine; out.pointing3DValid = 1; }
      }
    } catch (const std::exception&) {
      out = Result{};
      out.status = 4;
    }
    return out;
  }

private:
  KalmanVertexFitter fitter_;  // producer is edm::stream, so this instance is stream-local
};

#endif
