// pseudo vector class for SSE2

class v4df{
private:
	typedef double v2df __attribute__ ((vector_size(16)));
	typedef float  v4sf __attribute__ ((vector_size(16)));
	v2df l;
	v2df h;
public:
	v4df(){ h = l = (v2df){0.0, 0.0}; 
	}
	v4df(v2df _l, v2df _h) : l(_l), h(_h) {}
	v4df(double x) : l((v2df){x, x}), h((v2df){x, x}) {}
	v4df(double a, double b, double c, double d) :
		l((v2df){a, b}), h((v2df){c, d}) {}
	v4df(const double *p, bool aligned=true){
		if(aligned){
			l = *(v2df *)(p + 0);
			h = *(v2df *)(p + 2);
		}else{
			l = __builtin_ia32_loadupd(p + 0);
			h = __builtin_ia32_loadupd(p + 2);
		}
	}
	~v4df(){}

	friend v4df operator + (const v4df &, const v4df &);
	friend v4df operator - (const v4df &, const v4df &);
	friend v4df operator * (const v4df &, const v4df &);
	friend v4df operator / (const v4df &, const v4df &);
	friend v4df operator - (const v4df &v);

	v4df &operator += (const v4df &v){
		l += v.l;
		h += v.h;
		return *this;
	}
	v4df &operator -= (const v4df &v){
		l -= v.l;
		h -= v.h;
		return *this;
	}
	v4df &operator *= (const v4df &v){
		l *= v.l;
		h *= v.h;
		return *this;
	}
	v4df &operator /= (const v4df &v){
		l /= v.l;
		h /= v.h;
		return *this;
	}
	v4df inv(){
		return v4df((v2df){1.0, 1.0} / l, (v2df){1.0, 1.0} / h);
	}
	v4df sqrt(){
		return v4df(
			__builtin_ia32_sqrtpd(l),
			__builtin_ia32_sqrtpd(h));
	}
	double sum(){
		v2df tmp = l + h;
		return __builtin_ia32_vec_ext_v2df(tmp, 0)
			 + __builtin_ia32_vec_ext_v2df(tmp, 1);
	}

	template <int N>
	v4df & incl(){
		static const v2df val = {N, N};
		l += val;
		h += val;
		return *this;
	}
	v4df &operator ++(){
		return incl <1> ();
	}
	operator v4sf (){
		return __builtin_ia32_movlhps(
				__builtin_ia32_cvtpd2ps(l),
				__builtin_ia32_cvtpd2ps(h));
	}
};

inline v4df operator + (const v4df &u, const v4df &v){
	return v4df(u.l + v.l, u.h + v.h);
}
inline v4df operator - (const v4df &u, const v4df &v){
	return v4df(u.l - v.l, u.h - v.h);
}
inline v4df operator * (const v4df &u, const v4df &v){
	return v4df(u.l * v.l, u.h * v.h);
}
inline v4df operator / (const v4df &u, const v4df &v){
	return v4df(u.l / v.l, u.h / v.h);
}
inline v4df operator - (const v4df &v){
	return v4df(-v.l, -v.h);
}

