
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

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_SUPERPATH_H
#define PBRT_INTEGRATORS_SUPERPATH_H

// integrators/path.h*
#include "pbrt.h"
#include "integrator.h"
#include "lightdistrib.h"
#include <random>

namespace pbrt {

	// PathIntegrator Declarations
	class SuperPathIntegrator : public SamplerIntegrator {
	public:
		// PathIntegrator Public Methods
		SuperPathIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
			std::shared_ptr<Sampler> sampler,
			const Bounds2i &pixelBounds,
			const std::vector<std::shared_ptr<Film>>& _films,
			const std::shared_ptr<Material>& _pMaterial,
			Float rrThreshold = 1,
			const std::string &lightSampleStrategy = "spatial");

		void Preprocess(const Scene &scene, Sampler &sampler);
		Spectrum Li(const RayDifferential &ray, const Scene &scene,
			Sampler &sampler, MemoryArena &arena, int depth, const int64_t uCurSample = 0) const final;

	private:
		Spectrum UniformSampleOneLightEnv(const Interaction &it, const Scene &scene,
			MemoryArena &arena, Sampler &sampler,
			bool handleMedia = false,
			const Distribution1D *lightDistrib = nullptr,
			const int64_t uCurSample = 0) const;

		Spectrum EstimateDirectEnv(const Interaction &it, const Point2f &uShading,
			const Light &light, const Point2f &uLight,
			const Scene &scene, Sampler &sampler,
			MemoryArena &arena, bool handleMedia = false,
			bool specular = false) const;

		Vector3f ComputeWh(const Point2f &u);
		Float GetAlpha(const Point2f &u);

		Point2f Sample(std::unique_ptr<Sampler>& _pSamler, const Point2f &u);

	private:
		// PathIntegrator Private Data
		const int maxDepth;
		const Float rrThreshold;
		const std::string lightSampleStrategy;
		std::unique_ptr<LightDistribution> lightDistribution;
		std::shared_ptr<Material> m_pMaterial;
		MemoryArena m_Arena;
		SurfaceInteraction m_Isect;

		std::default_random_engine m_Generator;
		std::normal_distribution<float> m_NormalDist;
	};

	SuperPathIntegrator *CreateSuperPathIntegrator(const ParamSet &params,
		std::shared_ptr<Sampler> sampler,
		std::shared_ptr<const Camera> camera,
		const std::vector<std::shared_ptr<Film>>& _films,
		const std::shared_ptr<Material>& _pMaterial);

}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_SUPERPATH_H
