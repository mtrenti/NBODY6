// #include "v4df.h"

class v4sf {
public:
	typedef float  _v4sf __attribute__ ((vector_size(16)));
	typedef int    _v4si __attribute__ ((vector_size(16)));
	typedef double _v2df __attribute__ ((vector_size(16)));
private:
	_v4sf val;
public:
	v4sf(){
		val = (_v4sf){0.0f, 0.0f, 0.0f, 0.0f};
	}
	v4sf(_v4sf _val) : val(_val) {}
	v4sf(float f){
		val = (_v4sf){f, f, f, f};
	}
	v4sf(double f){
		val = (_v4sf){f, f, f, f};
	}
	v4sf(unsigned f){
		val = (_v4sf){f, f, f, f};
	}
	/*
	v4sf(const float &f, int){
		_v4sf tmp = __builtin_ia32_loadss(&f);
		val = __builtin_ia32_shufps(tmp, tmp, 0);
	}
	*/
	v4sf(float a, float b, float c, float d){
		val = (_v4sf){a, b, c, d};
	}
	v4sf(const float *p, bool aligned=1){
		if(aligned){
			val = *(_v4sf *)p;
		}else{
			val = __builtin_ia32_loadups(p);
		}
	}
	v4sf(const double *p, bool aligned=1){
		if(aligned){
			val = __builtin_ia32_movlhps(
					__builtin_ia32_cvtpd2ps(*(_v2df *)(p+0)),
					__builtin_ia32_cvtpd2ps(*(_v2df *)(p+2))
				  );
		}else{
			val = (_v4sf){p[0], p[1], p[2], p[3]};
		}
	}
	v4sf(_v4si ival) : val((_v4sf)ival) {}
	~v4sf() {}

	const v4sf &operator  = (v4sf a){
		val  = a.val;
		return *this;
	}
	const v4sf &operator += (v4sf a){
		val += a.val;
		return *this;
	}
	const v4sf &operator -= (v4sf a){
		val -= a.val;
		return *this;
	}
	const v4sf &operator *= (v4sf a){
		val *= a.val;
		return *this;
	}
	const v4sf &operator /= (v4sf a){
		val /= a.val;
		return *this;
	}
	const v4sf operator - () const{
		return v4sf(-val);
	}
	const v4sf operator + (v4sf a) const{
		return v4sf(val + a.val);
	}
	const v4sf operator - (v4sf a) const{
		return v4sf(val - a.val);
	}
	const v4sf operator * (v4sf a) const{
		return v4sf(val * a.val);
	}
	const v4sf operator / (v4sf a) const{
		return v4sf(val / a.val);
	}

	const v4sf operator != (v4sf a) const{
		return v4sf((_v4sf)__builtin_ia32_cmpneqps(val, a.val));
	}
	const v4sf operator == (v4sf a) const{
		return v4sf((_v4sf)__builtin_ia32_cmpeqps(val, a.val));
	}
	const v4sf operator & (v4sf a) const{
		return __builtin_ia32_andps(val, a.val);
	}

	// read-only for debugging
	float operator [] (int i){
		switch(i){
		case 0:
			return __builtin_ia32_vec_ext_v4sf(val, 0);
		case 1:
			return __builtin_ia32_vec_ext_v4sf(val, 1);
		case 2:
			return __builtin_ia32_vec_ext_v4sf(val, 2);
		case 3:
			return __builtin_ia32_vec_ext_v4sf(val, 3);
		default:
			return 0.0f;
		}
	}

	v4sf inv() const{
		return v4sf(1.0f)/(*this);
	}
	v4sf sqrt() const{
		return v4sf(__builtin_ia32_sqrtps(val));
	}
	v4sf rsqrt() const{
		v4sf x = val;
		v4sf y = __builtin_ia32_rsqrtps(x);
		return (v4sf(-0.5) * y) * (x*y*y + v4sf(-3.0));
	}

	operator v4df () const{
		_v2df vl = __builtin_ia32_cvtps2pd(val);
		_v2df vh = __builtin_ia32_cvtps2pd(
					__builtin_ia32_movhlps(val, val));
		return v4df(vl, vh);
	}

	operator _v4sf() const{
		return val;
	}
	_v2df as_v2df(){ // reinterpret cast
		return (_v2df)val;
	}
	_v4si as_v4si() const{
		return (_v4si)val;
	}

	_v2df dsum() const{
		_v2df vl = __builtin_ia32_cvtps2pd(val);
		_v2df vh = __builtin_ia32_cvtps2pd(
					__builtin_ia32_movhlps(val, val));
		return vl+vh;
	}
};
