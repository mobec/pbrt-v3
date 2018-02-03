
/*
pbrt source code is Copyright(c) 1998-2016
Matt Pharr, Greg Humphreys, and Wenzel Jakob.

This file is part of pbrt.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

// integrators/path.cpp*
#include "integrators/superpath.h"
#include "bssrdf.h"
#include "camera.h"
#include "film.h"
#include "interaction.h"
#include "paramset.h"
#include "scene.h"
#include "stats.h"
#include <fstream>
#include "samplers\random.h"

namespace pbrt {

	STAT_PERCENT("Integrator/Zero-radiance paths", zeroRadiancePaths, totalPaths);
	STAT_INT_DISTRIBUTION("Integrator/Path length", pathLength);

	// PathIntegrator Method Definitions
	SuperPathIntegrator::SuperPathIntegrator(int maxDepth,
		std::shared_ptr<const Camera> camera,
		std::shared_ptr<Sampler> sampler,
		const Bounds2i &pixelBounds,
		const std::vector<std::shared_ptr<Film>>& _films,
		const std::shared_ptr<Material>& _pMaterial,
		Float rrThreshold,
		const std::string &lightSampleStrategy)
		: SamplerIntegrator(camera, sampler, pixelBounds, _films),
		maxDepth(maxDepth),
		rrThreshold(rrThreshold),
		lightSampleStrategy(lightSampleStrategy),
		m_pMaterial(_pMaterial),
		m_NormalDist(0.f, 1.f)
	{
		m_Isect = {};
	}

	Vector3f SuperPathIntegrator::ComputeWh(const Point2f &u)
	{
		if (m_Isect.bsdf != nullptr)
		{
			return m_Isect.bsdf->GetWh(u);
		}

		return {};
	}

	Float SuperPathIntegrator::GetAlpha(const Point2f& u)
	{
		if (m_Isect.bsdf != nullptr)
		{
			return m_Isect.bsdf->GetAlpha(u);
		}

		return 1.f;
	}

	Point2f SuperPathIntegrator::Sample(std::unique_ptr<Sampler>& _pSamler, const Point2f &alpha)
	{
		float Norm = m_NormalDist(m_Generator);
		const Float fThreshold = fabsf(Norm); // deviation from normal

		const Vector3f vNormal = { 0, 0 , 1 };
		const auto goodsample = [&](const Point2f& u) 
		{
			Vector3f wh = ComputeWh(u);
			Float a = fabsf(Dot(wh, vNormal) - 1.f);
			return a  < fThreshold;
		};

		Point2f cur = _pSamler->Get2D();
		while (goodsample(cur) == false)
		{
			cur = _pSamler->Get2D();
		}
		return cur;
	}

	void SuperPathIntegrator::Preprocess(const Scene &scene, Sampler &sampler)
	{
		// is not really needed because we only have one light
		lightDistribution =	CreateLightSampleDistribution(lightSampleStrategy, scene);

		std::unique_ptr<Sampler> pLightSampler = sampler.Clone(0xB00B5); // 0xB00B5
		pLightSampler->StartPixel({});

		std::unique_ptr<Sampler> pScatterSampler = std::make_unique<RandomSampler>(sampler.samplesPerPixel, 0xA55A55);
		pScatterSampler->StartPixel({});

		uLights.resize(sampler.samplesPerPixel);
		uScatters.resize(sampler.samplesPerPixel);

		if (m_pMaterial != nullptr)
		{
			Bounds2i sampleBounds = camera->film->GetSampleBounds();
			Vector2i sampleExtent = sampleBounds.Diagonal();
			pbrt::Point2i pixel;
			pixel.x = sampleExtent.x / 2;
			pixel.y = sampleExtent.y / 2;

			// Initialize _CameraSample_ for current sample
			CameraSample cameraSample = pLightSampler->GetCameraSample(pixel);

			// Generate camera ray for current sample
			RayDifferential ray;
			Float rayWeight = camera->GenerateRayDifferential(cameraSample, &ray);
			ray.ScaleDifferentials(1 / std::sqrt((Float)pLightSampler->samplesPerPixel));
			bool foundIntersection = scene.Intersect(ray, &m_Isect);
			m_pMaterial->ComputeScatteringFunctions(&m_Isect, m_Arena, pbrt::TransportMode::Radiance, true);
		}

		for (size_t i = 0; i < sampler.samplesPerPixel; i++)
		{
			uLights[i] = pLightSampler->Get2D();
		}

		const Float alpha = GetAlpha({});
		uScatters[0] = { alpha + 1.f, 0.f}; // force wh == normal

		for (size_t i = 1; i < sampler.samplesPerPixel; i++)
		{
			if (m_Isect.bsdf != nullptr)
				uScatters[i] = Sample(pScatterSampler, uScatters[0]);
			else
				uScatters[i] = pScatterSampler->Get2D();
		}

		std::ofstream out("wh.txt");

		if (out.is_open())
		{			
			for (size_t i = 0; i <  sampler.samplesPerPixel; i++)
			{
				out << ComputeWh(uScatters[i]).z;					
				if(i+1 <sampler.samplesPerPixel) out << ",";
			}

			out.close();
		}

		//for (size_t i = 0; i < sampler.samplesPerPixel; i++)
		//{
		//	uScatters[i] = pScatterSampler->Get2D();
		//}
	}

	Spectrum SuperPathIntegrator::UniformSampleOneLightEnv(const Interaction &it, const Scene &scene,
		MemoryArena &arena, Sampler &sampler,
		bool handleMedia, const Distribution1D *lightDistrib, const int64_t uCurSample) const{
		ProfilePhase p(Prof::DirectLighting);
		// Randomly choose a single light to sample, _light_
		int nLights = int(scene.lights.size());
		if (nLights == 0) return Spectrum(0.f);
		int lightNum;
		Float lightPdf;
		if (lightDistrib) {
			lightNum = lightDistrib->SampleDiscrete(sampler.Get1D(), &lightPdf);
			if (lightPdf == 0) return Spectrum(0.f);
		}
		else {
			lightNum = std::min((int)(sampler.Get1D() * nLights), nLights - 1);
			lightPdf = Float(1) / nLights;
		}
		const std::shared_ptr<Light> &light = scene.lights[lightNum];
		Point2f uLight = sampler.Get2D();
		//Point2f uScattering = sampler.Get2D();
		//Point2f uLight = uLights[uCurSample];
		Point2f uScattering = uScatters[uCurSample];
		return EstimateDirectEnv(it, uScattering, *light, uLight,
			scene, sampler, arena, handleMedia) / lightPdf;
	}

	Spectrum SuperPathIntegrator::EstimateDirectEnv(const Interaction &it, const Point2f &uScattering,
		const Light &light, const Point2f &uLight,
		const Scene &scene, Sampler &sampler,
		MemoryArena &arena, bool handleMedia, bool specular) const {
		BxDFType bsdfFlags =
			specular ? BSDF_ALL : BxDFType(BSDF_ALL & ~BSDF_SPECULAR);
		Spectrum Ld(0.f);
		// Sample light source with multiple importance sampling
		Vector3f wi;
		Float lightPdf = 0, scatteringPdf = 0;
		VisibilityTester visibility;
		Spectrum Li = light.Sample_Li(it, uLight, &wi, &lightPdf, &visibility);

		VLOG(2) << "EstimateDirect uLight:" << uLight << " -> Li: " << Li << ", wi: "
			<< wi << ", pdf: " << lightPdf;
		if (lightPdf > 0 && !Li.IsBlack()) {
			// Compute BSDF or phase function's value for light sample
			Spectrum f;
			if (it.IsSurfaceInteraction()) {
				// Evaluate BSDF for light sampling strategy
				const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
				f = isect.bsdf->f(isect.wo, wi, bsdfFlags) *
					AbsDot(wi, isect.shading.n);
				scatteringPdf = isect.bsdf->Pdf(isect.wo, wi, bsdfFlags);
				VLOG(2) << "  surf f*dot :" << f << ", scatteringPdf: " << scatteringPdf;
			}
			else {
				// Evaluate phase function for light sampling strategy
				const MediumInteraction &mi = (const MediumInteraction &)it;
				Float p = mi.phase->p(mi.wo, wi);
				f = Spectrum(p);
				scatteringPdf = p;
				VLOG(2) << "  medium p: " << p;
			}
			if (!f.IsBlack()) {
				// Compute effect of visibility for light source sample
				if (handleMedia) {
					Li *= visibility.Tr(scene, sampler);
					VLOG(2) << "  after Tr, Li: " << Li;
				}
				else {
					if (!visibility.Unoccluded(scene)) {
						VLOG(2) << "  shadow ray blocked";
						Li = Spectrum(0.f);
					}
					else
						VLOG(2) << "  shadow ray unoccluded";
				}

				// Add light's contribution to reflected radiance
				//if (!Li.IsBlack()) {
				//	if (IsDeltaLight(light.flags))
				//		Ld += f * Li / lightPdf;
				//	else {
				//		Float weight =
				//			PowerHeuristic(1, lightPdf, 1, scatteringPdf);
				//		Ld += f * Li * weight / lightPdf;
				//	}
				//}
			}
		}

		// Sample BSDF with multiple importance sampling
		if (!IsDeltaLight(light.flags)) {
			Spectrum f;
			bool sampledSpecular = false;
			if (it.IsSurfaceInteraction()) {
				// Sample scattered direction for surface interactions
				BxDFType sampledType;
				const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
				f = isect.bsdf->Sample_f(isect.wo, &wi, uScattering, &scatteringPdf,
					bsdfFlags, &sampledType);
				f = AbsDot(wi, isect.shading.n);
				sampledSpecular = (sampledType & BSDF_SPECULAR) != 0;
			}
			else {
				// Sample scattered direction for medium interactions
				const MediumInteraction &mi = (const MediumInteraction &)it;
				Float p = mi.phase->Sample_p(mi.wo, &wi, uScattering);
				f = Spectrum(p);
				scatteringPdf = p;
			}
			VLOG(2) << "  BSDF / phase sampling f: " << f << ", scatteringPdf: " <<
				scatteringPdf;
			if (!f.IsBlack() && scatteringPdf > 0) {
				// Account for light contributions along sampled direction _wi_
				Float weight = 1;
				if (!sampledSpecular) {
					lightPdf = light.Pdf_Li(it, wi);
					if (lightPdf == 0) return Ld;
					weight = PowerHeuristic(1, scatteringPdf, 1, lightPdf);
				}

				// Find intersection and compute transmittance
				SurfaceInteraction lightIsect;
				Ray ray = it.SpawnRay(wi);
				Spectrum Tr(1.f);
				bool foundSurfaceInteraction =
					handleMedia ? scene.IntersectTr(ray, sampler, &lightIsect, &Tr)
					: scene.Intersect(ray, &lightIsect);

				// Add light contribution from material sampling
				Spectrum Li(0.f);
				if (foundSurfaceInteraction) {
					if (lightIsect.primitive->GetAreaLight() == &light)
						Li = lightIsect.Le(-wi);
				}
				else
					Li = light.Le(ray);

				if (!Li.IsBlack()) Ld += f * Li /** Tr * weight / scatteringPdf*/;
			}
		}
		return Ld;
	}

	Spectrum SuperPathIntegrator::Li(const RayDifferential &r, const Scene &scene,
		Sampler &sampler, MemoryArena &arena,
		int depth, const int64_t uCurSample) const {
		ProfilePhase p(Prof::SamplerIntegratorLi);
		Spectrum L(0.f), beta(1.f);
		RayDifferential ray(r);
		bool specularBounce = false;
		int bounces;
		// Added after book publication: etaScale tracks the accumulated effect
		// of radiance scaling due to rays passing through refractive
		// boundaries (see the derivation on p. 527 of the third edition). We
		// track this value in order to remove it from beta when we apply
		// Russian roulette; this is worthwhile, since it lets us sometimes
		// avoid terminating refracted rays that are about to be refracted back
		// out of a medium and thus have their beta value increased.
		Float etaScale = 1;

		for (bounces = 0;; ++bounces) {
			// Find next path vertex and accumulate contribution
			VLOG(2) << "Path tracer bounce " << bounces << ", current L = " << L
				<< ", beta = " << beta;

			// Intersect _ray_ with scene and store intersection in _isect_
			SurfaceInteraction isect;
			bool foundIntersection = scene.Intersect(ray, &isect);

			// Possibly add emitted light at intersection
			//if (bounces == 0 || specularBounce) {
			//	// Add emitted light at path vertex or from the environment
			//	if (foundIntersection) {
			//		L += beta * isect.Le(-ray.d);
			//		VLOG(2) << "Added Le -> L = " << L;
			//	}
			//	else {
			//		for (const auto &light : scene.infiniteLights)
			//			L += beta * light->Le(ray);
			//		VLOG(2) << "Added infinite area lights -> L = " << L;
			//	}
			//}

			if (bounces == 0 && foundIntersection == false)
			{
				Float rgb[] = { 0.f, 1.f, 0.f };
				L += beta * RGBSpectrum::FromRGB(rgb);
			}

			// Terminate path if ray escaped or _maxDepth_ was reached
			if (!foundIntersection || bounces >= maxDepth) break;

			// Compute scattering functions and skip over medium boundaries
			isect.ComputeScatteringFunctions(ray, arena, true);
			if (!isect.bsdf) {
				VLOG(2) << "Skipping intersection due to null bsdf";
				ray = isect.SpawnRay(ray.d);
				bounces--;
				continue;
			}

			const Distribution1D *distrib = lightDistribution->Lookup(isect.p);

			// Sample illumination from lights to find path contribution.
			// (But skip this for perfectly specular BSDFs.)
			if (isect.bsdf->NumComponents(BxDFType(BSDF_ALL & ~BSDF_SPECULAR)) >
				0) {
				++totalPaths;
				Spectrum Ld = beta * UniformSampleOneLightEnv(isect, scene, arena,
					sampler, false, distrib, uCurSample);
				VLOG(2) << "Sampled direct lighting Ld = " << Ld;
				if (Ld.IsBlack()) ++zeroRadiancePaths;
				CHECK_GE(Ld.y(), 0.f);
				L += Ld;
			}

			// Sample BSDF to get new path direction
			Vector3f wo = -ray.d, wi;
			Float pdf;
			BxDFType flags;
			Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf,
				BSDF_ALL, &flags);
			VLOG(2) << "Sampled BSDF, f = " << f << ", pdf = " << pdf;
			if (f.IsBlack() || pdf == 0.f) break;
			beta *= f * AbsDot(wi, isect.shading.n) / pdf;
			VLOG(2) << "Updated beta = " << beta;
			CHECK_GE(beta.y(), 0.f);
			DCHECK(!std::isinf(beta.y()));
			specularBounce = (flags & BSDF_SPECULAR) != 0;
			if ((flags & BSDF_SPECULAR) && (flags & BSDF_TRANSMISSION)) {
				Float eta = isect.bsdf->eta;
				// Update the term that tracks radiance scaling for refraction
				// depending on whether the ray is entering or leaving the
				// medium.
				etaScale *= (Dot(wo, isect.n) > 0) ? (eta * eta) : 1 / (eta * eta);
			}
			ray = isect.SpawnRay(wi);

			// Account for subsurface scattering, if applicable
			if (isect.bssrdf && (flags & BSDF_TRANSMISSION)) {
				// Importance sample the BSSRDF
				SurfaceInteraction pi;
				Spectrum S = isect.bssrdf->Sample_S(
					scene, sampler.Get1D(), sampler.Get2D(), arena, &pi, &pdf);
				DCHECK(!std::isinf(beta.y()));
				if (S.IsBlack() || pdf == 0) break;
				beta *= S / pdf;

				// Account for the direct subsurface scattering component
				L += beta * UniformSampleOneLight(pi, scene, arena, sampler, false,
					lightDistribution->Lookup(pi.p));

				// Account for the indirect subsurface scattering component
				Spectrum f = pi.bsdf->Sample_f(pi.wo, &wi, sampler.Get2D(), &pdf,
					BSDF_ALL, &flags);
				if (f.IsBlack() || pdf == 0) break;
				beta *= f * AbsDot(wi, pi.shading.n) / pdf;
				DCHECK(!std::isinf(beta.y()));
				specularBounce = (flags & BSDF_SPECULAR) != 0;
				ray = pi.SpawnRay(wi);
			}

			// Possibly terminate the path with Russian roulette.
			// Factor out radiance scaling due to refraction in rrBeta.
			Spectrum rrBeta = beta * etaScale;
			if (rrBeta.MaxComponentValue() < rrThreshold && bounces > 3) {
				Float q = std::max((Float).05, 1 - rrBeta.MaxComponentValue());
				if (sampler.Get1D() < q) break;
				beta /= 1 - q;
				DCHECK(!std::isinf(beta.y()));
			}
		}
		ReportValue(pathLength, bounces);
		return L;
	}

	SuperPathIntegrator *CreateSuperPathIntegrator(const ParamSet &params,
		std::shared_ptr<Sampler> sampler,
		std::shared_ptr<const Camera> camera,
		const std::vector<std::shared_ptr<Film>>& _films,
		const std::shared_ptr<Material>& _pMaterial) {
		int maxDepth = params.FindOneInt("maxdepth", 5);
		int np;
		const int *pb = params.FindInt("pixelbounds", &np);
		Bounds2i pixelBounds = camera->film->GetSampleBounds();
		if (pb) {
			if (np != 4)
				Error("Expected four values for \"pixelbounds\" parameter. Got %d.",
					np);
			else {
				pixelBounds = Intersect(pixelBounds,
					Bounds2i{ { pb[0], pb[2] },{ pb[1], pb[3] } });
				if (pixelBounds.Area() == 0)
					Error("Degenerate \"pixelbounds\" specified.");
			}
		}
		Float rrThreshold = params.FindOneFloat("rrthreshold", 1.);
		std::string lightStrategy =
			params.FindOneString("lightsamplestrategy", "spatial");
		return new SuperPathIntegrator(maxDepth, camera, sampler, pixelBounds,
			_films, _pMaterial, rrThreshold, lightStrategy);
	}

}  // namespace pbrt
