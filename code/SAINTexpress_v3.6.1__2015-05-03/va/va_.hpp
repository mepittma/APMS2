// begin(), end() for valarrays
namespace std{
	template<class T>
	inline T* begin(valarray<T>& a) {
		return &(a[0]);}

	template<class T>
	inline const T* begin(const valarray<T>& a) {
		return &(a[0]);}

	template<class T>
	inline T* end(valarray<T>& a) {
		return &(a[0]) + a.size();}

	template<class T>
	inline const T* end(const valarray<T>& a) {
		return &(a[0]) + a.size();}
}
