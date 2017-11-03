class space{};
class basis{};


class element_wise_addition : boost::function<array<double, 3>(array<double,3> const&, array<double,3> const&)>{
	operator(f
}
array<double, 3> element_wise_addition(array<double, 3> const& a, array<double,3> const& b){
	return (array<double, 3>){{a[0]+b[0], a[1]+b[1], a[2]+b[2]}};
}
array<double, 3> element_multiplication(array<double, 3> const& v, double const& s){
	return (array<double, 3>){{a[0]*s, a[1]*s, a[2]*s}};
}

template<class Set, class VectorAddition, class ScalarMultiplication>
class vector_space : private Set{
	static vector_space operator+(vector_space const& v1, vector_space const& v2){
		return VectorAddition(v1,v2);
	}
	static vector_space operator*(vector_space const& v, double const& s){
		return ScalarMultiplication(v,s);
	}
};

template<class Basis>
class covariant : vector_space<array<double, 3>, element_wise_addition, element_multiplication> {
};

template<class Basis>
class contravariant : array<double, 3>{
};

template<class Basis>
double operator*(contravariant<Basis> const& contra, covariant<Basis> const& co){
	return 
		contra[0]*co[0]+
		contra[1]*co[1]+
		contra[2]*co[2];
}

