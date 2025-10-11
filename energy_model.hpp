#pragma once

#include "miscs.hpp"

#define ENERGY_PARAM_NAMESPACE_BEGIN namespace energy { namespace turner2004 {
#define ENERGY_PARAM_NAMESPACE_END } }
#include "energy_param.hpp"
#undef ENERGY_PARAM_NAMESPACE_BEGIN
#undef ENERGY_PARAM_NAMESPACE_END

#define ENERGY_PARAM_NAMESPACE_BEGIN namespace energy { namespace turner1999 {
#define ENERGY_PARAM_NAMESPACE_END } }
#include "legacy_energy_param.hpp"
#undef ENERGY_PARAM_NAMESPACE_BEGIN
#undef ENERGY_PARAM_NAMESPACE_END

namespace energy {

enum class Model {
	Turner2004,
	Turner1999
};

struct Params{
	double temperature;
	double gas_constant;
	double k0;
	double kT;
	double lxc37;
	int ML_intern37;
	int ML_closing37;
	int ML_BASE37;
	int MAX_NINIO;
	int ninio37;
	int TerminalAU37;
	const int (*stack37)[NBPAIRS+1];
	const int *hairpin37;
	const int *bulge37;
	const int *internal_loop37;
	const int (*mismatchI37)[5][5];
	const int (*mismatch1nI37)[5][5];
	const int (*mismatch23I37)[5][5];
	const int (*mismatchH37)[5][5];
	const int (*mismatchM37)[5][5];
	const int (*mismatchExt37)[5][5];
	const int (*dangle5_37)[5];
	const int (*dangle3_37)[5];
	const int (*int11_37)[NBPAIRS+1][5][5];
	const int (*int21_37)[NBPAIRS+1][5][5][5];
	const int (*int22_37)[NBPAIRS+1][NBPAIRS+1][5][5][5][5];
	const char *Triloops;
	const int *Triloop37;
	const char *Tetraloops;
	const int *Tetraloop37;
	const char *Hexaloops;
	const int *Hexaloop37;
	bool has_special_hairpins;
	bool allow_mismatch_multi;
	bool allow_mismatch_external;
};

namespace detail{

// 実行時に作成する constexpr にしたため、かなり安全
constexpr Params make_turner2004(){
	return {
		.temperature = turner2004::temperature,
		.gas_constant = GASCONST,
		.k0 = K0,
		.kT = turner2004::kT,
		.lxc37 = turner2004::lxc37,
		.ML_intern37 = turner2004::ML_intern37,
		.ML_closing37 = turner2004::ML_closing37,
		.ML_BASE37 = turner2004::ML_BASE37,
		.MAX_NINIO = turner2004::MAX_NINIO,
		.ninio37 = turner2004::ninio37,
		.TerminalAU37 = turner2004::TerminalAU37,
		.stack37 = turner2004::stack37,
		.hairpin37 = turner2004::hairpin37,
		.bulge37 = turner2004::bulge37,
		.internal_loop37 = turner2004::internal_loop37,
		.mismatchI37 = turner2004::mismatchI37,
		.mismatch1nI37 = turner2004::mismatch1nI37,
		.mismatch23I37 = turner2004::mismatch23I37,
		.mismatchH37 = turner2004::mismatchH37,
		.mismatchM37 = turner2004::mismatchM37,
		.mismatchExt37 = turner2004::mismatchExt37,
		.dangle5_37 = turner2004::dangle5_37,
		.dangle3_37 = turner2004::dangle3_37,
		.int11_37 = turner2004::int11_37,
		.int21_37 = turner2004::int21_37,
		.int22_37 = turner2004::int22_37,
		.Triloops = turner2004::Triloops,
		.Triloop37 = turner2004::Triloop37,
		.Tetraloops = turner2004::Tetraloops,
		.Tetraloop37 = turner2004::Tetraloop37,
		.Hexaloops = turner2004::Hexaloops,
		.Hexaloop37 = turner2004::Hexaloop37,
		.has_special_hairpins = true,
		.allow_mismatch_multi = true,
		.allow_mismatch_external = true
	};
}

constexpr Params make_turner1999(){
	return {
		.temperature = turner1999::temperature,
		.gas_constant = GASCONST,
		.k0 = K0,
		.kT = turner1999::kT,
		.lxc37 = turner1999::lxc37,
		.ML_intern37 = turner1999::ML_intern37,
		.ML_closing37 = turner1999::ML_closing37,
		.ML_BASE37 = turner1999::ML_BASE37,
		.MAX_NINIO = turner1999::MAX_NINIO,
		.ninio37 = turner1999::ninio37,
		.TerminalAU37 = turner1999::TerminalAU37,
		.stack37 = turner1999::stack37,
		.hairpin37 = turner1999::hairpin37,
		.bulge37 = turner1999::bulge37,
		.internal_loop37 = turner1999::internal_loop37,
		.mismatchI37 = turner1999::mismatchI37,
		.mismatch1nI37 = turner1999::mismatch1nI37,
		.mismatch23I37 = turner1999::mismatch23I37,
		.mismatchH37 = turner1999::mismatchH37,
		.mismatchM37 = nullptr,
		.mismatchExt37 = nullptr,
		.dangle5_37 = turner1999::dangle5_37,
		.dangle3_37 = turner1999::dangle3_37,
		.int11_37 = turner1999::int11_37,
		.int21_37 = turner1999::int21_37,
		.int22_37 = turner1999::int22_37,
		.Triloops = nullptr,
		.Triloop37 = nullptr,
		.Tetraloops = nullptr,
		.Tetraloop37 = nullptr,
		.Hexaloops = nullptr,
		.Hexaloop37 = nullptr,
		.has_special_hairpins = false,
		.allow_mismatch_multi = false,
		.allow_mismatch_external = false
	};
}

} // namespace detail

inline const Params& get_params(Model model){
	static const Params turner2004_params = detail::make_turner2004();
	static const Params turner1999_params = detail::make_turner1999();
	switch(model){
	case Model::Turner1999:
		return turner1999_params;
	case Model::Turner2004:
	default:
		return turner2004_params;
	}
}

} // namespace energy
