#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -std=c++14 -O3 -Wall -Wextra -fmax-errors=2 `#-Wfatal-errors` -D_TEST_MPI3_COMMUNICATOR $0x.cpp -o $0x.x && time mpirun -np 1 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
//  (C) Copyright Alfredo A. Correa 2018.
#ifndef MPI3_COMMUNICATOR_HPP
#define MPI3_COMMUNICATOR_HPP

#include "../mpi3/communication_mode.hpp"
#include "../mpi3/equality.hpp"
#include "../mpi3/generalized_request.hpp"
#include "../mpi3/handle.hpp"
#include "../mpi3/info.hpp"
#include "../mpi3/message.hpp"
#include "../mpi3/operation.hpp"
#include "../mpi3/port.hpp"
#include "../mpi3/request.hpp"
#include "../mpi3/status.hpp"
#include "../mpi3/type.hpp"
#include "../mpi3/error.hpp"

#include "../mpi3/detail/basic_communicator.hpp"
#include "../mpi3/detail/buffer.hpp"
#include "../mpi3/detail/datatype.hpp"
#include "../mpi3/detail/iterator.hpp"
#include "../mpi3/detail/value_traits.hpp"
#include "../mpi3/detail/strided.hpp"
#include "../mpi3/detail/package.hpp"

#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

#include <boost/optional.hpp> // TODO replace by std::optional (c++17)
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_array.hpp>

#define BOOST_PACKAGE_ARCHIVE_SOURCE

#include <boost/archive/detail/common_oarchive.hpp>
#include <boost/archive/detail/common_iarchive.hpp>
#include <boost/archive/basic_streambuf_locale_saver.hpp>
#include <boost/archive/archive_exception.hpp>
#include <boost/archive/detail/auto_link_archive.hpp>
#include <boost/archive/detail/abi_prefix.hpp> // must be the last header

#include <boost/serialization/item_version_type.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/is_bitwise_serializable.hpp>

#include <boost/mpl/placeholders.hpp>

// use this to avoid need for linking -lserialization
#ifdef _MAKE_BOOST_SERIALIZATION_HEADER_ONLY
#if BOOST_VERSION < 106600 and BOOST_VERSION > 106000
#include "../mpi3/serialization_hack/singleton.cpp"
#endif
#if BOOST_VERSION < 105900
#define BOOST_ARCHIVE_DECL
#define BOOST_SERIALIZATION_DECL
#endif
#include "../mpi3/serialization_hack/archive_exception.cpp"
#include "../mpi3/serialization_hack/extended_type_info.cpp"
#include "../mpi3/serialization_hack/extended_type_info_typeid.cpp"
#include "../mpi3/serialization_hack/basic_archive.cpp"
#include "../mpi3/serialization_hack/basic_iarchive.cpp"
#include "../mpi3/serialization_hack/basic_oarchive.cpp"
#include "../mpi3/serialization_hack/basic_iserializer.cpp"
#include "../mpi3/serialization_hack/basic_oserializer.cpp"
#endif

#include "../mpi3/package_archive.hpp"

#include<cassert>
#include<iostream>
#include<iterator> // iterator_traits
#include<limits>
#include<map>
#include<numeric> // std::accumulate
#include<string>
#include<thread>
#include<type_traits>
#include<vector>
#include<type_traits> // is_same
		

namespace boost{
namespace mpi3{

#if !defined(OPEN_MPI) || (OMPI_MAJOR_VERSION < 2)
#define OMPI_COMM_TYPE_NODE     MPI_COMM_TYPE_SHARED
#define OMPI_COMM_TYPE_HWTHREAD MPI_COMM_TYPE_SHARED
#define OMPI_COMM_TYPE_CORE     MPI_COMM_TYPE_SHARED
#define OMPI_COMM_TYPE_L1CACHE  MPI_COMM_TYPE_SHARED
#define OMPI_COMM_TYPE_L2CACHE  MPI_COMM_TYPE_SHARED
#define OMPI_COMM_TYPE_L3CACHE  MPI_COMM_TYPE_SHARED
#define OMPI_COMM_TYPE_SOCKET   MPI_COMM_TYPE_SHARED
#define OMPI_COMM_TYPE_NUMA     MPI_COMM_TYPE_SHARED
#define OMPI_COMM_TYPE_BOARD    MPI_COMM_TYPE_SHARED
#define OMPI_COMM_TYPE_HOST     MPI_COMM_TYPE_SHARED
#define OMPI_COMM_TYPE_CU       MPI_COMM_TYPE_SHARED
#define OMPI_COMM_TYPE_CLUSTER  MPI_COMM_TYPE_SHARED
#endif

enum class communicator_type : int{
	node = OMPI_COMM_TYPE_NODE,
	hw_thread = OMPI_COMM_TYPE_HWTHREAD,
	core = OMPI_COMM_TYPE_CORE,
	l1_cache = OMPI_COMM_TYPE_L1CACHE,
	l2_cache = OMPI_COMM_TYPE_L2CACHE,
	l3_cache = OMPI_COMM_TYPE_L3CACHE,
	socket = OMPI_COMM_TYPE_SOCKET,
	numa = OMPI_COMM_TYPE_NUMA,
	board = OMPI_COMM_TYPE_BOARD,
	host = OMPI_COMM_TYPE_HOST,
	cu = OMPI_COMM_TYPE_CU,
	cpu = OMPI_COMM_TYPE_CU,
	cluster = OMPI_COMM_TYPE_CLUSTER 
};

template<int N = 10> struct overload_priority : overload_priority<N-1>{
//	using overload_priority<N-1>::overload_priority;
};
template<> struct overload_priority<0>{};

class environment;
class group;

using boost::optional;

template<class T = void>
struct window;

struct error_handler;
struct keyval;

template<class T>
struct shm_pointer;

template<class T>
struct pointer;

using address = MPI_Aint;
using intptr_t = MPI_Aint;
using size_t = MPI_Aint;

struct request;

struct send_request;
struct receive_request;


struct FILE;

struct process;

struct ostream;

struct message_header{
	int tag;
	
};

struct graph_communicator;
struct shared_communicator; // intracommunicator

class communicator : protected detail::basic_communicator{
	friend struct detail::package;
protected:
	bool is_null() const{return MPI_COMM_NULL == impl_;}
	friend class mpi3::environment;
	equality compare(communicator const& other) const{
		equality ret = boost::mpi3::unequal;
		MPI_Comm_compare(impl_, other.impl_, reinterpret_cast<int*>(&ret));
		return ret;
	}
public:

//	using detail::basic_communicator::send;
//	using detail::basic_communicator::send_n;
	using detail::basic_communicator::send_receive_n;
	using detail::basic_communicator::matched_probe;
	template<class It, typename Size>
	auto send_n(
		It first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size count,
		int dest, int tag
	){
		int s = MPI_Send(
			detail::data(first), count, 
			detail::basic_datatype<typename std::iterator_traits<It>::value_type>{},
			dest, tag, impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot send");
	}
	template<class It, typename Size>
	auto isend_n(
		It first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size count,
		int dest, int tag
	){
		mpi3::request r;
		int s = MPI_Isend(
			detail::data(first), count, 
			detail::basic_datatype<typename std::iterator_traits<It>::value_type>{},
			dest, tag, impl_, &r.impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot send");
		return r;
	}
	template<class It, typename Size>
	void send_n(
		It first, 
			detail::forward_iterator_tag,
			detail::value_unspecified_tag,
		Size count,
		int dest, int tag
	){
		detail::package p(*this);
		package_oarchive poa(p);
		while(count--) poa << *first++;
		send_n(p.begin(), p.size(), dest, tag); //	p.send(dest, tag);
	}
	template<class It, typename Size>
	auto isend_n(It first, Size count, int dest, int tag = 0){
		return isend_n(
			first, 
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			count,
			dest, tag
		);
	}
	template<class It, typename Size>
	auto isend_n(
		It first, 
			detail::forward_iterator_tag,
			detail::value_unspecified_tag,
		Size count, 
		int dest, int tag
	){
		detail::package p(*this);
		package_oarchive poa(p);
		while(count--) poa << *first++;
		return isend_n(p.begin(), p.size(), dest, tag); //	p.send(dest, tag);
	}
	template<class It, typename Size>
	void send_n(It first, Size count, int dest, int tag = 0){
		return send_n(
			first, 
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			count,
			dest, tag
		);
	}
	template<class It>
	auto send(
		It first, It last,
			detail::random_access_iterator_tag, 
			detail::value_unspecified_tag,
		int dest, int tag
	){
		return send_n(first, std::distance(first, last), dest, tag);
	}
	template<class It>
	auto send(
		It first, It last,
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		int dest, int tag
	){
		return send_n(first, std::distance(first, last), dest, tag);
	}
	template<class It>
	auto send(
		It first, It last,
			detail::input_iterator_tag, 
			detail::basic_tag,
		int dest, int tag
	){
		mpi3::vector<typename std::iterator_traits<It>::value_type> buffer(first, last);
		return send_n(buffer.begin(), buffer.size(), dest, tag);
	}
	template<class It>
	auto send(
		It first, It last,
			detail::input_iterator_tag,
			detail::value_unspecified_tag,
		int dest, int tag
	){
		detail::package p(*this);
		package_oarchive poa(p);
		while(first!=last) poa << *first++;
		send_n(p.begin(), p.size(), dest, tag); //	p.send(dest, tag);
	}
	template<class It>
	auto isend(
		It first, It last,
			detail::random_access_iterator_tag, 
			detail::value_unspecified_tag,
		int dest, int tag
	){
		return isend_n(first, std::distance(first, last), dest, tag);
	}
	template<class It>
	auto send(It first, It last, int dest, int tag = 0){
		return send(
			first, last,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			dest, tag
		);
	}
	template<class It>
	auto isend(It first, It last, int dest, int tag = 0){
		return isend(
			first, last,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			dest, tag
		);
	}
	using detail::basic_communicator::basic_communicator;
	communicator(communicator const&) = default;
	communicator(communicator&&) = default;
	communicator() = default;
	communicator& operator=(communicator const& other){
		communicator tmp(other);
		swap(tmp);
		return *this;
	}
	communicator& operator=(communicator&& other){
		communicator tmp(std::move(other));
		swap(tmp);
		return *this;
	}
	bool operator==(communicator const& other) const{
		auto eq = compare(other);
		return (eq == equality::identical) or (eq == equality::congruent);
	}
	bool operator!=(communicator const& other) const{return not (*this == other);}
	explicit operator bool() const{return not is_null();}
	impl_t operator&() const{return impl_;}
	~communicator(){
		if(impl_ != MPI_COMM_WORLD and impl_ != MPI_COMM_NULL and impl_ != MPI_COMM_SELF){
			MPI_Comm_disconnect(&impl_); //	or MPI_Comm_free(&impl_);?
		}
	}
	int size() const{
		if(is_null()) return 0;//throw std::runtime_error("size called on null communicator");
		int size = -1; 
		int s = MPI_Comm_size(impl_, &size);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot get size"); 
		return size;
	}
	bool empty() const{return size() == 0;}
	void abort(int errorcode = 0) const{MPI_Abort(impl_, errorcode);}
	bool is_intercommunicator() const{
		int flag = false;
		MPI_Comm_test_inter(impl_, &flag);
		return flag;
	}
	communicator split(int color, int key) const{
		communicator ret;
		int s = MPI_Comm_split(impl_, color, key, &ret.impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot split communicator");
		if(ret) ret.name(name() + std::to_string(color));// + std::to_string(key));
		return ret;
	}
	communicator split(int color = MPI_UNDEFINED) const{return split(color, rank());}
	shared_communicator split_shared(int key = 0) const;
	shared_communicator split_shared(communicator_type t, int key = 0) const;
	int remote_size() const{
		int ret = -1;
		int s = MPI_Comm_remote_size(impl_, &ret);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot remote size");
		return ret;
	}
	communicator reversed() const{
		return split(0, size() - rank());
	}
	int cartesian_map(std::vector<int> const& dims, std::vector<int> const& periods) const{
		assert( dims.size() == periods.size() );
		int ret;
		int s = MPI_Cart_map(impl_, dims.size(), dims.data(), periods.data(), &ret);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot map"};
		return ret;
	}
	int cartesian_map(std::vector<int> const& dimensions) const{
		return cartesian_map(dimensions, std::vector<int>(dimensions.size(), 0));
	}
	template<class T>
	window<T> make_window(mpi3::size_t n); // Win_allocate
	template<class T = void>
	window<T> make_window(T* base, mpi3::size_t n);
	template<class T = void>
	window<T> make_window();

	pointer<void> malloc(MPI_Aint size) const;
	template<class T = void>
	void deallocate_shared(pointer<T> p);
	template<class T = void>
	void deallocate(pointer<T>& p, MPI_Aint size = 0);
	void free(pointer<void>& p) const;

	bool similar(communicator const& other) const{
		return compare(other) == equality::similar;
	}
	template<class Vector>//, typename = typename std::enable_if<std::is_same<decltype(Vector{}.data()), int*>{}>::type>
	communicator subcomm(Vector const& v) const{
		MPI_Group old_g;
		MPI_Comm_group(impl_, &old_g);
		MPI_Group new_g;
		MPI_Group_incl(old_g, v.size(), v.data(), &new_g);
		communicator ret;
		MPI_Comm_create(impl_, new_g, &ret.impl_);
		return ret;
	}
	communicator subcomm(std::initializer_list<int> l) const{
		return subcomm(std::vector<int>(l));
	}

	enum class topology{undefined = MPI_UNDEFINED, graph = MPI_GRAPH, cartesian = MPI_CART};
//	topology topo() const{return static_cast<topology>(call<&MPI_Topo_test>());}
	int rank() const{
		int rank = -1;
		if(impl_ == MPI_COMM_NULL) throw std::runtime_error("rank on null communicator");
		int s = MPI_Comm_rank(impl_, &rank);
		if(s != MPI_SUCCESS) MPI_Comm_call_errhandler(impl_, s);
		return rank;
	}
	int right() const{return (rank() + 1) % size();}
	int left() const{
		int left = rank() - 1;
		if(left < 0) left = size() - 1;
		return left;
	}
	int next(int n = 1) const{
		assert(rank() + n < size());
		return rank() + n;
	}
	int prev(int n = 1) const{
		assert(rank() - n > 0);
		return rank() - n;
	}
	communicator accept(port const& p, int root = 0) const{
		communicator ret;
		MPI_Comm_accept(p.name_.c_str(), MPI_INFO_NULL, root, impl_, &ret.impl_);
		return ret;
	}
	void barrier() const{MPI_Barrier(impl_);}
	communicator connect(port const& p, int root = 0) const{
		communicator ret;
		MPI_Comm_connect(p.name_.c_str(), MPI_INFO_NULL, root, impl_, &ret.impl_);
		return ret;
	}
	bool root() const{return rank() == 0;}
	void set_error_handler(error_handler const& eh);
	error_handler get_error_handler() const;
	template<typename T>
	void set_attribute(keyval const& kv, int idx, T const& val);
	void delete_attribute(keyval const& kv, int idx);
	template<typename T>
	void get_attribute(keyval const& kv, int idx, T& val);
	template<typename T>
	T const& get_attribute_as(keyval const& kv, int idx);
	bool has_attribute(keyval const& kv, int idx);

	process operator[](int i);

	template<typename T>
	T const& attribute_as(int TAG) const{
		int flag = 0;
		T* p = nullptr;
		int status = MPI_Comm_get_attr(impl_, TAG, &p, &flag);
		if(status != MPI_SUCCESS) throw std::runtime_error{"cannot get attr"};
		assert(flag);
		return *p;
	}

	void call_error_handler(int errorcode){
		int status = MPI_Comm_call_errhandler(impl_, errorcode);
		if(status != MPI_SUCCESS) throw std::runtime_error{"cannot call error handler"};
	}
	void error(mpi3::error const& e){
		int s = MPI_Comm_call_errhandler(impl_, static_cast<int>(e));
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot call error handler"};
	}


//	mpi3::group group() const;
	
	communicator operator/(int n) const{
		int color = rank()*n/size();
		return split(color);
	}
	communicator operator%(int n) const{
		int color = rank()%n;
		return split(color);
	}
	communicator operator<(int n) const{return split((rank() < n)?0:MPI_UNDEFINED);}
	communicator operator<=(int n) const{return split((rank() <= n)?0:MPI_UNDEFINED);}
	communicator operator>(int n) const{return split((rank() > n)?0:MPI_UNDEFINED);}
	communicator operator>=(int n) const{return split((rank() >= n)?0:MPI_UNDEFINED);}
	communicator operator==(int n) const{return split((rank() == n)?0:MPI_UNDEFINED);}

	template<class T>
	void send_value(T const& t, int dest, int tag = 0){
		send(std::addressof(t), std::addressof(t) + 1, dest, tag);
	}
	template<class T>
	auto isend_value(T const& t, int dest, int tag = 0){
		return isend(std::addressof(t), std::addressof(t) + 1, dest, tag);
	}

	template<class T, std::size_t N>
	void send_value(T(&t)[N], int dest, int tag = 0){
		send(std::addressof(t[0]), std::addressof(t[N-1]) + 1, dest, tag);
	}

//	template<class T>
//	void receive_value(T& t, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
//		return receive(std::addressof(t), std::addressof(t) + 1, source, tag);
//	}
//	template<class T>
//	auto ireceive_value(T& t, int source = MPI_ANY_TAG, int tag = MPI_ANY_TAG){
//		return ireceive(std::addressof(t), std::addressof(t) + 1, source, tag);
//	}
//	template<class T, std::size_t N>
//	void receive_value(T(&t)[N], int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
//		receive(std::addressof(t[0]), std::addressof(t[N]), source, tag);
//	}
#if 0
	template<
		class ContiguousIterator, 
		typename Size, 
		class datatype = detail::basic_datatype<typename std::iterator_traits<ContiguousIterator>::value_type>
	//	, typename = std::enable_if_t<detail::iterator_stride<ContiguousIterator>==1>
	> 
	void send_n_aux(std::random_access_iterator_tag, ContiguousIterator it, Size count, int dest, int tag = 0) const{
		MPI_Send(std::addressof(*it), count, datatype{}, dest, tag, impl_);
	}
	template<
		class RandomAccessIterator, 
		class Size, class value_type = typename std::iterator_traits<RandomAccessIterator>::value_type, 
		class datatype = detail::basic_datatype<value_type>
	>
	void send_n_aux(std::random_access_iterator_tag, RandomAccessIterator I, Size count, int dest, int tag = 0) const{
		MPI_Send(std::addressof(*I), count, datatype{}, dest, tag, impl_);
	}
	template<class InputIterator, class Size>
	void send_n_aux(std::input_iterator_tag, InputIterator I, Size count, int dest, int tag = 0) const{
		for(Size i = 0; i != count; ++i) send_value(*I++, dest, tag);
	}
	template<class InputIterator, class... A, class Category = typename std::iterator_traits<InputIterator>::iterator_category>
	void send_n(InputIterator I, A&&... a) const{send_n_aux(Category(), I, std::forward<A>(a)...);}
#endif

	template<class ContIt, class Size>
//	send_request 
	request send_init_n(ContIt first, Size count, int dest, int tag = 0);
	template<class ContIt>
	request
//	send_request 
	send_init(ContIt first, ContIt last, int dest, int tag = 0);

	template<class ContIt, class Size>
	receive_request receive_init_n(ContIt first, Size count, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG);
	template<class ContIt>
	receive_request receive_init(ContIt first, ContIt last, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG);

	template<typename InputIt, typename Size, 
		class value_type = typename std::iterator_traits<InputIt>::value_type, 
		class datatype = detail::basic_datatype<value_type>
	>
	auto pack_n_aux(std::random_access_iterator_tag, InputIt first, Size count, char* out, int max_size
	){
		assert(0);
		int position = 0;
		MPI_Pack(detail::data(first), count, datatype{}, out, max_size, &position, impl_);
		return out + position;
	}

	template<class InputIt, class Size, typename V = typename std::iterator_traits<InputIt>::value_type, class datatype = detail::basic_datatype<V>>
	auto unpack_from_n(InputIt first, Size count, char const* buffer){
		int position = 0;
		using detail::data;
	//	assert( count % pack_size<V>() == 0 );
		int s = MPI_Unpack(buffer, count*pack_size<V>(), &position, data(first), count, datatype{}, impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot packaga unpack");
		return buffer + position;
	}
//	template<class InputIt, typename Size>
//	auto unpack_n(InputIt first, Size count, char const* buffer){
//		return unpack_from_n(first, count, buffer);
//	}
//	template<class V>
//	auto unpack_value(V& v, char const* buffer){
//		return unpack_n(&v, 1, buffer);
//	}
//	template<class V>
//	auto pack_value(V const& v, char const* buffer){
//		return pack_n(&v, 1, buffer);
//	}
	template<class It>
	auto pack(
		It first, It last, 
			detail::random_access_iterator_tag,
		uvector<detail::packed>& b, int pos
	){
		return pack_n(first, std::distance(first, last), b, pos);
	}
	template<class It>
	auto pack(It first, It last, uvector<detail::packed>& b, int pos){
		return pack(
			first, last,
				detail::iterator_category_t<It>{},
			b, pos
		);
	}
	template<class It, class Size>
	auto send_receive_replace_n(
		It first, Size size, 
		int dest, int source, // = MPI_ANY_SOURCE, 
		int sendtag = 0, int recvtag = MPI_ANY_TAG
	){
		using value_type = typename std::iterator_traits<It>::value_type;
		return send_receive_replace_n(
			first, 
				detail::iterator_category_t<It>{}, 
				detail::value_category_t<value_type>{},
			size,
			dest, source, sendtag, recvtag
		);
	}
	template<class It, typename Size>
	auto send_receive_replace_n(
		It first, 
			detail::random_access_iterator_tag,
			detail::basic_tag,
		Size count, int dest, int source, int sendtag, int recvtag
	){
		status ret;
		int s = MPI_Sendrecv_replace(
			detail::data(first), count, 
			detail::basic_datatype<typename std::iterator_traits<It>::value_type>{}, 
			dest, sendtag, source, recvtag, impl_, &ret.impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot sendrecv_replace"};
		return ret;
	}
	template<class It, typename Size>
	auto send_receive_n(
		It first, Size count, int dest, 
		It d_first, Size d_count, int source, 
		int sendtag = 0, int recvtag = MPI_ANY_TAG
	){
		return send_receive_n(
			first, count, dest, 
			d_first, d_count, source,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{}, 
			sendtag, recvtag
		);
	}
	template<class It, typename Size>
	auto send_receive_n(
		It first, Size count, int dest, 
		It d_first, Size d_count, int source,
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		int sendtag, int recvtag
	){
		status ret;
		int s = MPI_Sendrecv(
			detail::data(first), count, 
			detail::basic_datatype<typename std::iterator_traits<It>::value_type>{},
			dest, sendtag, 
			detail::data(d_first), d_count,
			detail::basic_datatype<typename std::iterator_traits<It>::value_type>{},
			source, recvtag,
			impl_,
			&ret.impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot sendrecv"};
		return ret;
	}
	template<class It, typename Size, typename... Meta>
	auto send_receive_replace_n(
		It first, 
			detail::forward_iterator_tag,
			detail::value_unspecified_tag,
		Size count, 
		int dest, int source,
		int sendtag, int recvtag
	){
		detail::package p(*this);
		package_oarchive poa(p);
		auto first_copy = first;
		while(count--) poa << *first_copy++;
		auto s = p.size();
		send_receive_replace_n(&s, 1, dest, source, sendtag, recvtag);
		detail::package p2(*this);
		p2.resize(s);
		send_receive_n(
			p.begin(), p.size(), dest, 
			p2.begin(), p2.size(), 
			source, sendtag, recvtag
		);
		package_iarchive pia(p2);
		while(p2) pia >> *first++;
		return first;
	}
	template<class It, typename Size>
	auto send_receive_replace_n(
		It first, 
			detail::forward_iterator_tag,
			detail::basic_tag,
		Size count, int dest, int source, int sendtag, int recvtag
	){
		uvector<typename std::iterator_traits<It>::value_type> v(count);
		std::copy_n(first, count, v.begin());
		send_receive_replace_n(v.begin(), v.size(), dest, source, sendtag, recvtag);
		return std::copy_n(v.begin(), v.size(), first);
	}
	template<class It, class Size>
	auto send_receive_n(
		It first, Size size, 
		int dest, int source, // = MPI_ANY_SOURCE, 
		int sendtag = 0, int recvtag = MPI_ANY_TAG
	){
		return send_receive_replace_n(
			first, size,
			dest, source, sendtag, recvtag
		);
	}
public:
	template<class It>
	auto send_receive(
		It first, It last, int dest, 
		It d_first, It d_last, int source, 
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		int sendtag, int recvtag
	){
		return send_receive_n(
			first, std::distance(first, last), dest,
			d_first, std::distance(d_first, d_last), source,
			sendtag, recvtag
		);
	}
	template<class It>
	auto send_receive(
		It first, It last, int dest, 
		It d_first, It d_last, int source, 
		int sendtag = 0, int recvtag = MPI_ANY_TAG
	){
		return send_receive(
			first, last, dest,
			d_first, d_last, source,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{}, 
			sendtag, recvtag
		);
	}
/*
		int size = comm_.pack_size<V>()*count;
		auto old_buffer_size = buffer_.size();
		buffer_.resize(old_buffer_size + size);
		auto curr = buffer_.data() + old_buffer_size;
		int position = -1;
		int outcount = this->pack_size<V>()*count;
		using detail::data;
//		auto end = comm_.pack_n(data(first), count, curr);
		MPI_Pack(detail::data(first), count, datatype{}, out, outcount, &position, impl_);
		end = out + position;
		assert(end == std::addressof(*buffer_.end()));
		in_ = buffer_.size();
	//	in_ += end - curr;
	}
*/

public:
	status probe(int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		status s;
		MPI_Probe(source, tag, impl_, &s.impl_);
		return s;
	}
	optional<status> iprobe(int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		status s;
		int flag;
		MPI_Iprobe(source, tag, impl_, &flag, &s.impl_);
		if(not flag) return optional<status>();
		return optional<status>(s);
	}
	template<class It>
	auto send_receive_replace(
		It first, It last, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		int dest, int source, int sendtag, int recvtag
	){
		return send_receive_replace_n(
			first, std::distance(first, last), 
			dest, source, sendtag, recvtag
		);		
	}
	template<class It>
	auto send_receive_replace(
		It first, It last, 
			detail::forward_iterator_tag,
			detail::basic_tag,
		int dest, int source, int sendtag, int recvtag
	){
		mpi3::vector<typename std::iterator_traits<It>::value_type> buffer(first, last);
		send_receive_replace_n(buffer.begin(), buffer.size(), dest, source, sendtag, recvtag);
		return std::copy(buffer.begin(), buffer.end(), first);
	}
	template<class It>
	auto send_receive_replace(
		It first, It last, 
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		int dest, int source, int sendtag, int recvtag
	){
		return send_receive_replace_n(
			first, std::distance(first, last), 
			dest, source, sendtag, recvtag
		);
	}
	template<class It>
	auto send_receive_replace(It first, It last, int dest, int source = MPI_ANY_SOURCE, int sendtag = 0, int recvtag = MPI_ANY_TAG){
		return send_receive_replace(
			first, last,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			dest, source, sendtag, recvtag
		);
	}
	template<class It>
	auto send_receive(It first, It last, int dest, int source = MPI_ANY_SOURCE, int sendtag = 0, int recvtag = MPI_ANY_TAG){
		return send_receive_replace(first, last, dest, source, sendtag, recvtag);
	}
/*	template<class It, class... Args>
	auto send_receive_replace_category(std::random_access_iterator_tag const&, 
		It first, It last, Args&&... args){
		return send_receive_replace_n(first, std::distance(first, last), std::forward<Args>(args)...);
	}*/
/*
	template<class T1, class T2>
	auto send_receive_value(
		T1 const& t1, int dest, 
		T2&       t2, int source = MPI_ANY_SOURCE
	){
		return send_receive_n(&t1, 1, dest, &t2, 1, source);
	}
	template<
		class IOIterator, class Size, class... A, 
		class CategoryIO = typename std::iterator_traits<IOIterator>::iterator_category,
		class datatype = detail::basic_datatype<typename std::iterator_traits<IOIterator>::value_type>
	>
	auto send_receive_n_aux(
		std::random_access_iterator_tag, IOIterator IO, Size count, int dest, int source, 
		int sendtag = 0, int recvtag = MPI_ANY_SOURCE
	) const{
		status ret;
		MPI_Sendrecv_replace(std::addressof(*IO), count, datatype{}, dest, sendtag, source, recvtag, impl_, &ret.impl_);
		return ret;
	}
*/
/*	template<
		class IOIterator, class Size, class... A, 
		class CategoryIO = typename std::iterator_traits<IOIterator>::iterator_category
	>
	auto send_receive_n(IOIterator IO, Size count, int dest, int source, A&&... a) const{
		return send_receive_n_aux(CategoryIO(), IO, count, dest, source, std::forward<A>(a)...);
	}*/

//	template<class T>
//	auto send_receive_value(T& t, int dest, int source = MPI_ANY_SOURCE, int sendtag = 0, int recvtag = MPI_ANY_TAG) const{
//		return send_receive_n(&t, 1, dest, source, sendtag, recvtag);
//	}

	void send_packed_n(void const* begin, int n, int dest, int tag = 0){
		std::cout << "sending packet of size " << n << std::endl;
		MPI_Send(const_cast<void*>(begin), n, MPI_PACKED, dest, tag, impl_);
	}
//	void send_value(package const& p, int dest, int tag = 0);
//	template<class Char = char>
//	auto send_packed(Char const* begin, Char const* end, int dest, int tag = 0){
//		return send_packed_n(begin, (char*)end - (char*)begin, dest, tag);
//	}
	auto receive_packed_n(void* begin, int n, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG) const{
		status ret;
		MPI_Recv(begin, n, MPI_PACKED, source, tag, impl_, &ret.impl_);
		return ret;
	}
	auto receive_packed(void* begin, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		MPI_Status status;
		MPI_Message msg;
		int count = -1;
		MPI_Mprobe(source, tag, impl_, &msg, &status);
		MPI_Get_count(&status, MPI_PACKED, &count);
		MPI_Mrecv(begin, count, MPI_PACKED, &msg, MPI_STATUS_IGNORE);
	//	auto n = probe(source, tag).count<char>();
	//	receive_packed_n(begin, n, source, tag);
		return (void*)((char*)begin + count);
	}
	template<class It, typename Size>
	auto receive_n(
		It dest, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size count, 
		int source, int tag
	){
		status sta;
		int s = MPI_Recv(
			detail::data(dest), count, 
			detail::basic_datatype<typename std::iterator_traits<It>::value_type>{},
			source, tag, impl_, &sta.impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("receive_n");
		return dest + count;
	}
	template<class It, typename Size>
	auto ireceive_n(
		It dest, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size count, 
		int source, int tag
	){
		mpi3::request r;
		int s = MPI_Irecv(
			detail::data(dest), count, 
			detail::basic_datatype<typename std::iterator_traits<It>::value_type>{},
			source, tag, impl_, &r.impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("receive_n");
		return r;
	//	return dest + count;
	}
	template<class It, typename Size>
	auto receive_n(
		It dest, 
			detail::forward_iterator_tag,
			detail::value_unspecified_tag,
		Size count, 
		int source, int tag
	){
		detail::package p(*this);
		p.receive(source, tag);
		package_iarchive pia(p);
		while(count--) pia >> *dest++;
		return dest;
	}
/*	template<class It, typename Size>
	auto ireceive_n(
		It dest, 
			detail::forward_iterator_tag,
			detail::value_unspecified_tag,
		Size count, 
		int source, int tag
	){
		detail::package p(*this);
		p.receive(source, tag);
		package_iarchive pia(p);
		while(count--) pia >> *dest++;
		return dest;
	}*/
/*	template<class It, typename Size>
	receive_request ireceive_n(
		It dest, 
			detail::forward_iterator_tag,
			detail::value_unspecified_tag,
		Size count, 
		int source, int tag
	){
		detail::package p(*this);
		p.receive(source, tag);
		package_iarchive pia(p);
		while(count--) pia >> *dest++;
		return dest;
	}*/
	template<class It, typename Size>
	auto receive_n(It dest, Size n, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		return receive_n(
			dest, 
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			n,
			source, tag
		);
	}
	template<class It, typename Size>
	mpi3::request ireceive_n(
		It dest, Size n, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG
	){
		return ireceive_n(
			dest,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			n,
			source, tag
		);
	}
	template<class It>
	auto receive(
		It dest, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		int source, int tag
	){
		match m = matched_probe(source, tag);
		auto count = m.count<typename std::iterator_traits<It>::value_type>();
		m.receive_n(dest, count);
		return dest + count;
	}
	template<class It>
	auto receive(
		It dest, 
			detail::forward_iterator_tag,
			detail::value_unspecified_tag,
		int source, int tag
	){
		detail::package p(*this);
		p.receive(source, tag);
		package_iarchive pia(p);
		while(p) pia >> *dest++;
		return dest;
	}
	template<class It>
	auto receive(
		It dest, 
			detail::forward_iterator_tag,
			detail::basic_tag,
		int source, int tag
	){
	//	match m = matched_probe(source, tag);
	//	auto count = m.count<typename std::iterator_traits<It>::value_type>();
	//	mpi3::uvector<typename std::iterator_traits<It>::value_type> v(count);
		return matched_probe(source, tag).receive_n(dest);//, count);
	}
	template<class It>
	auto receive(It dest, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		return receive(
			dest,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			source, tag
		);
	}
	template<class It>
	auto receive(
		It d_first, It d_last, 
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		int source, int tag
	){
		return receive_n(d_first, std::distance(d_first, d_last), source, tag);
	}
	template<class It>
	auto receive(
		It d_first, It d_last, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		int source, int tag
	){
		return receive_n(std::addressof(*d_first), std::distance(d_first, d_last), source, tag);
	//	return std::copy(buffer.begin(), buffer.end(), d_first);
	}

	template<class It>
	auto receive(
		It d_first, It d_last, 
			detail::forward_iterator_tag,
			detail::basic_tag,
		int source, int tag
	){
		mpi3::uvector<typename std::iterator_traits<It>::value_type> buffer(std::distance(d_first, d_last));
		receive_n(buffer.begin(), buffer.size(), source, tag);
		return std::copy(buffer.begin(), buffer.end(), d_first);
	}
	class ir_req{
		boost::mpi3::status query(){
			boost::mpi3::status ret;
			ret.set_source(MPI_UNDEFINED);
			ret.set_tag(MPI_UNDEFINED);
			ret.set_cancelled();
			ret.set_elements<char>(0);
			return ret;
		}
		void free(){
			std::cout << "free" << std::endl;
		}
		void cancel(int complete){
			std::cout << "cancel " << complete << std::endl;
		}
	};
	template<class It>
	struct receive_args{ 
		communicator* commP; 
		It d_first;
	//	It d_last;
		int source;
		int tag; 
		MPI_Request* requestP; 
	};
	struct receive_state{
		int cancelled = 0;
		int source = MPI_UNDEFINED;
		int tag = MPI_UNDEFINED;
	};
	template<class It>
	inline static void* receive_thread(void* ptr){
		receive_args<It>* args = (receive_args<It>*)ptr;
		args->commP->receive(args->d_first, args->source, args->tag);//, /*args->d_last,*/ );
		MPI_Grequest_complete(*args->requestP);
		::free(ptr);
		return NULL;
	}
	inline static int query_fn(void* extra_state, MPI_Status *status){
		receive_state* rs = (receive_state*)extra_state;
		/* always send just one int */ 
		MPI_Status_set_elements(status, MPI_INT, 1);
		/* can never cancel so always true */ 
		MPI_Status_set_cancelled(status, rs->cancelled); 
		/* choose not to return a value for this */
		status->MPI_SOURCE = rs->source; 
		/* tag has not meaning for this generalized request */ 
		status->MPI_TAG = rs->tag; 
		/* this generalized request never fails */ 
		return MPI_SUCCESS; 
	}
	inline static int free_fn(void* extra_state){ 
		/* this generalized request does not need to do any freeing */ 
		/* as a result it never fails here */
		::free(extra_state);
		return MPI_SUCCESS; 
	}
	inline static int cancel_fn(void* /*extra_state*/, int complete) 
	{ 
		/* This generalized request does not support cancelling. 
		   Abort if not already done.  If done then treat as if cancel failed. */ 
		if(!complete){
			fprintf(stderr, "Cannot cancel generalized request - aborting program\n"); 
			MPI_Abort(MPI_COMM_WORLD, 99); 
		} 
		return MPI_SUCCESS;
	}
	template<class It>
	auto ireceive(It d_first, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		// based on http://liinwww.ira.uka.de/courses/spprakt/mpi2-html-doc/node157.html
		mpi3::request ret; /*	receive_args<It>* args = (receive_args<It>*)::malloc(sizeof(receive_args<It>)); args->commP = this; args->d_first = d_first; //	args->d_last = d_last; args->source = source; args->tag = tag; args->requestP = &ret.impl_;*/
		receive_state* rs = (receive_state*)::malloc(sizeof(receive_state));
		rs->cancelled = 0;
		rs->source = source;
		rs->tag = tag;
		MPI_Grequest_start(query_fn, free_fn, cancel_fn, rs, &ret.impl_);//args->requestP);
		std::thread( //	static_cast<void*(*)(void*)>(receive_thread<It>), args
			[this, d_first, source, tag, &ret](){
				this->receive(d_first, source, tag); //	receive_args<It>* args = (receive_args<It>*)ptr; //	args->commP->receive(args->d_first, args->source, args->tag);//, /*args->d_last,*/ );
				MPI_Grequest_complete(ret.impl_); //	MPI_Grequest_complete(*args->requestP); //	::free(ptr);
			}
		).detach();	//	t.detach(); //	pthread_t thread; //	pthread_create(&thread, NULL, static_cast<void*(*)(void*)>(receive_thread<It>), args); //	pthread_detach(thread);
		return ret;		
	}
	template<class It>
	auto ireceive(
		It d_first, It d_last, 
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		int source, int tag
	){
		return ireceive_n(d_first, std::distance(d_first, d_last), source, tag);
	}
	template<class It>
	auto receive(It d_first, It d_last, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		return receive(
			d_first, d_last,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			source, tag
		);
	}
	template<typename InputIt, typename Size, class V = typename std::iterator_traits<InputIt>::value_type>
	auto bsend_init_n(
		InputIt first, Size count, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		int dest, int tag
	){
		request ret;
		int s = MPI_Bsend_init(
			std::addressof(*first), count, detail::basic_datatype<V>{},
			dest, tag, impl_, &ret.impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot bsend init"};
		return ret;
	}
	template<typename It, typename Size>
	auto bsend_init_n(It first, Size count, int dest, int tag = 0){
		return bsend_init_n(
			first, count,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			dest, tag
		);				
	}
	template<class It>
	auto bsend_init(
		It first, It last, 
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		int dest, int tag
	){
		return bsend_init_n(first, std::distance(first, last), dest, tag);
	}
	template<class It>
	auto bsend_init(It first, It last, int dest, int tag = 0){
		bsend_init(
			first, last,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			dest, tag
		);				
	}
	template<class InputIterator, class category = typename std::iterator_traits<InputIterator>::iterator_category>
	auto bsend(InputIterator It1, InputIterator It2, int dest, int tag = 0){
		return send(buffered_communication_mode{}, blocking_mode{}, It1, It2, dest, tag);
	}
	template<class InputIt, class V = typename std::iterator_traits<InputIt>::value_type>
	auto dynamic_receive(InputIt first, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
	//	auto count = probe(source, tag).count<V>();
	//	return receive(first, first + count, source, tag);
		MPI_Status status;
	    MPI_Message msg;
    	int count = -1;
    	MPI_Mprobe(source, tag, impl_, &msg, &status);
    	MPI_Get_count(&status, detail::basic_datatype<V>{}, &count);
    //	int* buffer = (int*)malloc(sizeof(int) * count);
    	using detail::data;
    	MPI_Mrecv(data(first), count, detail::basic_datatype<V>{}, &msg, MPI_STATUS_IGNORE);
	}

/*	template<class InputIterator>
	auto receive(InputIterator first, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		return dynamic_receive(first, source, tag);
	}*/ 

	template<class Iterator, class category = typename std::iterator_traits<Iterator>::iterator_category>
	auto breceive(Iterator It1, Iterator It2, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		return receive(buffered_communication_mode{}, blocking_mode{}, It1, It2, source, tag);
	}
	template<class InputIterator, class category = typename std::iterator_traits<InputIterator>::iterator_category>
	auto ibsend(InputIterator It1, InputIterator It2, int dest, int tag = 0){
		return send(buffered_communication_mode{}, nonblocking_mode{}, It1, It2, dest, tag);
	}
	template<class Iterator, class category = typename std::iterator_traits<Iterator>::iterator_category>
	auto ibreceive(Iterator It1, Iterator It2, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		return receive(buffered_communication_mode{}, nonblocking_mode{}, It1, It2, source, tag);
	}

	template<class InputIterator, class category = typename std::iterator_traits<InputIterator>::iterator_category>
	auto ssend(InputIterator It1, InputIterator It2, int dest, int tag = 0){
		return send(synchronous_communication_mode{}, blocking_mode{}, It1, It2, dest, tag);
	}
	template<class Iterator, class category = typename std::iterator_traits<Iterator>::iterator_category>
	auto sreceive(Iterator It1, Iterator It2, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		return receive(synchronous_communication_mode{}, blocking_mode{}, It1, It2, source, tag);
	}
	template<class InputIterator, class category = typename std::iterator_traits<InputIterator>::iterator_category>
	auto issend(InputIterator It1, InputIterator It2, int dest, int tag = 0){
		return send(synchronous_communication_mode{}, nonblocking_mode{}, It1, It2, dest, tag);
	}
	template<class CommunicationMode, class It1, class Size>
	auto isend_n(
		CommunicationMode, 
		It1 first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size count, int dest, int tag
	){
		mpi3::request ret;
		CommunicationMode{}.ISend(
			std::addressof(*first), count, detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{}, dest, tag, impl_, &ret.impl_
		);
		return ret;
	}
	template<class It1, class Size>
	auto issend_n(It1 first, Size count, int dest, int tag = 0){
		return isend_n(
			synchronous_communication_mode{}, 
			first, 
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			count, 
			dest, tag
		);
	}
	template<class Iterator, class category = typename std::iterator_traits<Iterator>::iterator_category>
	auto isreceive(Iterator It1, Iterator It2, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		return receive(synchronous_communication_mode{}, nonblocking_mode{}, It1, It2, source, tag);
	}

	template<class InputIterator, class category = typename std::iterator_traits<InputIterator>::iterator_category>
	auto rsend(InputIterator It1, InputIterator It2, int dest, int tag = 0){
		return send(ready_communication_mode{}, blocking_mode{}, It1, It2, dest, tag);
	}
	template<class Iterator, class category = typename std::iterator_traits<Iterator>::iterator_category>
	auto rreceive(Iterator It1, Iterator It2, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		return receive(ready_communication_mode{}, blocking_mode{}, It1, It2, source, tag);
	}
	template<class InputIterator, class category = typename std::iterator_traits<InputIterator>::iterator_category>
	auto irsend(InputIterator It1, InputIterator It2, int dest, int tag = 0){
		return send(ready_communication_mode{}, nonblocking_mode{}, It1, It2, dest, tag);
	}
	template<class Iterator, class category = typename std::iterator_traits<Iterator>::iterator_category>
	auto irreceive(Iterator It1, Iterator It2, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG){
		return receive(ready_communication_mode{}, nonblocking_mode{}, It1, It2, source, tag);
	}
	template<class CommunicationMode, class BlockingMode, class InputIterator>
	auto send_category(CommunicationMode cm, BlockingMode bm, std::input_iterator_tag, 
		InputIterator first, InputIterator last, int dest, int tag
	);
	template<class CommunicationMode, class BlockingMode, class InputIterator, class category = typename std::iterator_traits<InputIterator>::iterator_category>
	auto send(CommunicationMode cm, BlockingMode bm, InputIterator It1, InputIterator It2, int dest, int tag = 0){
		return send_category(cm, bm, category{}, It1, It2, dest, tag);
	}
private:
/*	{
		auto it = mpi3::type::registered.find(typeid(V)); assert(it != boost::mpi3::type::registered.end());
		int s = cm.Send(detail::data(first), n, it->second.impl_, dest, tag, impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot send random access iterators");
	}*/
//	{
//		auto p = make_package();
//		assert(0);
#if 0
		auto it = boost::mpi3::type::registered.find(typeid(ValueType));
		if(it != boost::mpi3::type::registered.end()){
			int s = cm.Send(detail::data(first), n, it->second.impl_, dest, tag, impl_);
			if(s != MPI_SUCCESS) throw std::runtime_error("cannot send random access iterators");
		}else if(std::is_pod<ValueType>{}){
#ifdef BOOST_MPI3_DISALLOW_AUTOMATIC_POD_COMMUNICATION
			throw std::logic_error("cannot send unregistered type " + std::string(typeid(ValueType).name()) + " (pod comm disallowed)");
#endif
			// allow seamless receive of pod types
			int s = cm.Send(detail::data(first), n*sizeof(ValueType), MPI_BYTE, dest, tag, impl_);
			if(s != MPI_SUCCESS) throw std::runtime_error("cannot send random access iterators");
		}else throw std::logic_error("cannot send unregistered non pod type " + std::string(typeid(ValueType).name()));
#endif
//	}

//	template<class CommunicationMode, class ContiguousIterator, class Size, class V = typename std::iterator_traits<ContiguousIterator>::value_type>
//	void receive_n_randomaccess_contiguous_builtin(CommunicationMode cm, blocking_mode, std::false_type, ContiguousIterator first, Size n, int dest, int tag);
/*	{
		status stts;
		auto it = mpi3::type::registered.find(typeid(V)); assert(it != boost::mpi3::type::registered.end());
		int s = cm.Recv(detail::data(first), n, it->second.impl_, dest, tag, impl_, &stts.impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot send random access iterators");
	}*/
#if 0
{
		auto it = boost::mpi3::type::registered.find(typeid(ValueType));
		if(it != boost::mpi3::type::registered.end()){
			int s = cm.Recv(detail::data(first), n, it->second.impl_, dest, tag, impl_, MPI_STATUS_IGNORE);
			if(s != MPI_SUCCESS) throw std::runtime_error("cannot send random access iterators");
		}else if(std::is_pod<ValueType>{}){
#ifdef BOOST_MPI3_DISALLOW_AUTOMATIC_POD_COMMUNICATION
			throw std::logic_error("cannot receive not registered type " + std::string(typeid(ValueType).name()) + " (pod comm disallowed)");
#endif
			// allow seamless receive of pod types
			int s = cm.Recv(detail::data(first), sizeof(ValueType)*n, MPI_BYTE, dest, tag, impl_, MPI_STATUS_IGNORE);
			if(s != MPI_SUCCESS) throw std::runtime_error("cannot send random access iterators");
		}else throw std::logic_error("cannot receive not registered non pod type " + std::string(typeid(ValueType).name()));
	}
#endif

	template<class CommunicationMode, class ContiguousIterator, class Size, class ValueType = typename std::iterator_traits<ContiguousIterator>::value_type, class datatype = typename detail::basic_datatype<ValueType> >
	request send_n_randomaccess_contiguous_builtin(CommunicationMode cm, nonblocking_mode, std::true_type, ContiguousIterator first, Size n, int dest, int tag){
		request r;
		int s = cm.ISend(detail::data(first), n, datatype{}, dest, tag, impl_, &r.impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot send random access iterators"};
		return r;
	}
	template<class CommunicationMode, class ContiguousIterator, class Size, class ValueType = typename std::iterator_traits<ContiguousIterator>::value_type, class datatype = typename detail::basic_datatype<ValueType> >
	request receive_n_randomaccess_contiguous_builtin(CommunicationMode cm, nonblocking_mode, std::true_type, ContiguousIterator first, Size n, int dest, int tag){
		request r;
		int s = cm.IRecv(detail::data(first), n, datatype{}, dest, tag, impl_, &r.impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot send random access iterators");
		return r;
	}

#if 0
public:
	template<class It1, class It2>
	auto all_gather(It1 first, It1 last, It2 d_first){
		return all_gather_builtinQ(
			std::integral_constant<bool, 
				detail::is_builtin_type<typename std::iterator_traits<It1>::value_type>{} and
				detail::is_builtin_type<typename std::iterator_traits<It2>::value_type>{}
			>{}, first, last, d_first
		);
	}
private:
	template<class It1, class It2>
	auto all_gather_builtinQ(std::true_type, It1 first, It1 last, It2 d_first){
		return all_gather_builtin(
			detail::iterator_category<It1>{}, 
			detail::iterator_category<It2>{}, first, last, d_first
		);
	}
	template<class It1, class It2, 
		class sendtype = detail::basic_datatype<typename std::iterator_traits<It1>::value_type>,
		class recvtype = detail::basic_datatype<typename std::iterator_traits<It2>::value_type>
	>
	auto all_gather_builtin(detail::contiguous_iterator_tag, detail::contiguous_iterator_tag, It1 first, It1 last, It2 d_first){
		int s = MPI_Allgather( detail::data(first), std::distance(first, last), 
			sendtype::value,
			detail::data(d_first), std::distance(first, last),
			recvtype::value,
			impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot scatter");
	}
	template<class Iterator1, class Iterator2>
	auto all_gather_builtinQ(std::false_type, Iterator1 first, Iterator2 last, Iterator1 d_first)
//	{ TODO implement }
	;
#endif

#if 0
public:
	template<class InputIt, class OutputIt, class InputCategory = typename std::iterator_traits<OutputIt>::iterator_category, class OutputCategory = typename std::iterator_traits<OutputIt>::iterator_category >
	auto all_gather(InputIt first, InputIt last, OutputIt d_first, OutputIt d_last){
		return all_gather_category(InputCategory{}, OutputCategory{}, first, last, d_first, d_last);
	}
private:
	template<class InputIt, class OutputIt>
	auto all_gather_category(std::random_access_iterator_tag, std::random_access_iterator_tag, InputIt first, InputIt last, OutputIt d_first, OutputIt d_last){
		return all_gather_n_randomaccess(first, std::distance(first, last), d_first, std::distance(d_first, d_last));
	}
	template<class RandomAccessIt1, class RandomAccessIt2, class Size>
	auto all_gather_n_randomaccess(RandomAccessIt1 first, Size count, RandomAccessIt2 d_first, Size d_count){
		return all_gather_n_randomaccess_contiguity(detail::is_contiguous<RandomAccessIt1>{}, detail::is_contiguous<RandomAccessIt2>{}, first, count, d_first, d_count);
	}
	template<class ContiguousIt1, class ContiguousIt2, class Size>
	auto all_gather_n_randomaccess_contiguity(std::true_type, std::true_type, ContiguousIt1 first, Size count, ContiguousIt2 d_first, Size d_count){
		return all_gather_n_randomaccess_contiguous_builtin(
			detail::is_builtin_type<typename std::iterator_traits<ContiguousIt1>::value_type>{}, 
			detail::is_builtin_type<typename std::iterator_traits<ContiguousIt2>::value_type>{}, 
			first, count, d_first, d_count
		);
	}
	template<class ContiguousIt1, class ContiguousIt2, class Size, class ValueType1 = typename std::iterator_traits<ContiguousIt1>::value_type, class ValueType2 = typename std::iterator_traits<ContiguousIt2>::value_type, class datatype1 = typename detail::basic_datatype<ValueType1>, class datatype2 = typename detail::basic_datatype<ValueType2> >
	auto all_gather_n_randomaccess_contiguous_builtin(std::true_type, std::true_type, ContiguousIt1 first, Size count, ContiguousIt2 d_first, Size d_count){
		int s = MPI_Allgather(detail::data(first), count, datatype1{}, detail::data(d_first), d_count, datatype2{}, impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot all gather");
	}
#endif

private:
	template<class It1, typename Size, class It2>
	auto all_to_all_n(
		It1 first,
			detail::contiguous_iterator_tag,
			detail::basic_tag, 
		Size count, 
		It2 d_first,
			detail::contiguous_iterator_tag,
			detail::basic_tag
	){
		int s = MPI_Alltoall( 
			detail::data(first), count, 
			detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{},
			detail::data(d_first), count,
			detail::basic_datatype<typename std::iterator_traits<It2>::value_type>{},
			impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot all2all"};
		return d_first + count;
	}
public:
	template<class It1, typename Size, class It2>
	auto all_to_all_n(It1 first, Size count, It2 d_first){
		return 
			all_to_all_n(
				first, 
					detail::iterator_category_t<It1>{},
					detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
				count,
				d_first,
					detail::iterator_category_t<It2>{},
					detail::value_category_t<typename std::iterator_traits<It2>::value_type>{}
			);
	}
public:
private:
/*	template<class BlockingMode, class InputIt, class OutputIt, class InputCategory = typename std::iterator_traits<OutputIt>::iterator_category, class OutputCategory = typename std::iterator_traits<OutputIt>::iterator_category >
	auto gather(BlockingMode bm, InputIt first, InputIt last, OutputIt d_first, OutputIt d_last, int root = 0){
		return gather_category(bm, InputCategory{}, OutputCategory{}, first, last, d_first, d_last, root);
	}
	template<class BlockingMode, class InputIt, class OutputIt>
	auto gather_category(BlockingMode bm, std::random_access_iterator_tag, std::random_access_iterator_tag, InputIt first, InputIt last, OutputIt d_first, OutputIt d_last, int root){
		return gather_n_randomaccess(bm, first, std::distance(first, last), d_first, std::distance(d_first, d_last), root);
	}
	template<class BlockingMode, class RandomAccessIt1, class RandomAccessIt2, class Size>
	auto gather_n_randomaccess(BlockingMode bm, RandomAccessIt1 first, Size count, RandomAccessIt2 d_first, Size d_count, int root){
		return gather_n_randomaccess_contiguity(detail::is_contiguous<RandomAccessIt1>{}, detail::is_contiguous<RandomAccessIt2>{}, first, count, d_first, d_count, root);
	}
	template<class BlockingMode, class ContiguousIt1, class ContiguousIt2, class Size>
	auto gather_n_randomaccess_contiguity(BlockingMode bm, std::true_type, std::true_type, ContiguousIt1 first, Size count, ContiguousIt2 d_first, Size d_count, int root){
		return gather_n_randomaccess_contiguous_builtin(bm,
			detail::is_basic<typename std::iterator_traits<ContiguousIt1>::value_type>{}, 
			detail::is_basic<typename std::iterator_traits<ContiguousIt2>::value_type>{}, 
			first, count, d_first, d_count, root
		);
	}
	template<class BlockingMode, class ContiguousIt1, class ContiguousIt2, class Size>
	auto gather_n_randomaccess_contiguity(BlockingMode bm, std::false_type, std::true_type, ContiguousIt1 first, Size count, ContiguousIt2 d_first, Size d_count, int root){
		return gather_n_randomaccess_noncontiguous_contiguous_strided(bm, 
			detail::is_strided<ContiguousIt1>{}, detail::is_strided<ContiguousIt2>{},
			first, count, d_first, d_count, root
		);
	}
	template<class BlockingMode, class ContiguousIt1, class ContiguousIt2, class Size>
	auto gather_n_randomaccess_noncontiguous_contiguous_strided(BlockingMode bm, std::true_type, std::false_type, ContiguousIt1 first, Size count, ContiguousIt2 d_first, Size d_count, int root){
		return gather_n_randomaccess_noncontiguous_contiguous_strided_nonstrided_builtin(bm, 
			detail::is_basic<typename std::iterator_traits<ContiguousIt1>::value_type>{}, 
			detail::is_basic<typename std::iterator_traits<ContiguousIt2>::value_type>{}, 
			first, count, d_first, d_count, root
		);
	}
	template<class BlockingMode, class ContiguousIt1, class ContiguousIt2, class Size, class ValueType1 = typename std::iterator_traits<ContiguousIt1>::value_type, class ValueType2 = typename std::iterator_traits<ContiguousIt2>::value_type, class datatype1 = typename detail::basic_datatype<ValueType1>, class datatype2 = typename detail::basic_datatype<ValueType2> >
	auto gather_n_randomaccess_noncontiguous_contiguous_strided_nonstrided_builtin(blocking_mode,std::true_type, std::true_type, ContiguousIt1 first, Size count, ContiguousIt2 d_first, Size d_count, int root){
		mpi3::type datatype1_strided = mpi3::type(datatype1{}).vector(count, 1, first.stride());
		int s = MPI_Gather(detail::data(first.base()), 1, datatype1{}, detail::data(d_first), d_count, datatype2{}, root, impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot all gather");
	}
	template<class BlockingMode, class ContiguousIt1, class ContiguousIt2, class Size, class ValueType1 = typename std::iterator_traits<ContiguousIt1>::value_type, class ValueType2 = typename std::iterator_traits<ContiguousIt2>::value_type, class datatype1 = typename detail::basic_datatype<ValueType1>, class datatype2 = typename detail::basic_datatype<ValueType2> >
	auto gather_n_randomaccess_noncontiguous_contiguous_strided_nonstrided_builtin(nonblocking_mode,std::true_type, std::true_type, ContiguousIt1 first, Size count, ContiguousIt2 d_first, Size d_count, int root){
		mpi3::type datatype1_strided = mpi3::type(datatype1{}).vector(count, 1, first.stride());
		request ret;
		int s = MPI_Igather(detail::data(first.base()), 1, datatype1{}, detail::data(d_first), d_count, datatype2{}, root, impl_, &ret.impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot Igather");
		return ret;
	}
*/
/*	template<class ContiguousIt1, class ContiguousIt2, class Size, class ValueType1 = typename std::iterator_traits<ContiguousIt1>::value_type, class ValueType2 = typename std::iterator_traits<ContiguousIt2>::value_type, class datatype1 = typename detail::basic_datatype<ValueType1>, class datatype2 = typename detail::basic_datatype<ValueType2> >
	auto gather_n_randomaccess_contiguous_builtin(std::true_type, std::true_type, ContiguousIt1 first, Size count, ContiguousIt2 d_first, Size d_count, int root){
		int s = MPI_Gather(detail::data(first), count, datatype1{}, detail::data(d_first), d_count, datatype2{}, root, impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot gather");
	}*/
private:
	template<class It1, class It2>
	auto scatter_builtinQ(std::true_type, It1 first, It1 last, It2 d_first, int root){
		return scatter_builtin(
			detail::iterator_category<It1>{}, 
			detail::iterator_category<It2>{}, first, last, d_first, root
		);
	}
	template<class It1, class It2, 
		class sendtype = detail::basic_datatype<typename std::iterator_traits<It1>::value_type>,
		class recvtype = detail::basic_datatype<typename std::iterator_traits<It2>::value_type>
	>
	auto scatter_builtin(detail::contiguous_iterator_tag, detail::contiguous_iterator_tag, It1 first, It1 last, It2 d_first, int root){
		int s = MPI_Scatter( detail::data(first), std::distance(first, last), 
			sendtype::value,
			detail::data(d_first), std::distance(first, last),
			recvtype::value,
			root,
			impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot scatter");
	}
	template<class Iterator1, class Iterator2>
	auto scatter_builtinQ(std::false_type, Iterator1 first, Iterator2 last, Iterator1 d_first, int root)
//	{ TODO implement }
	;
public:
	template<class T>
	auto broadcast_value(T& t, int root = 0){
		return broadcast_n(std::addressof(t), 1, root);
	}
#if 0
public:
	template<class T>
	auto broadcast_value(T& t, int root = 0){return broadcast_n(std::addressof(t), 1, root);}
	template<class InputIt>
	auto broadcast(InputIt first, InputIt last, int root = 0){
		return broadcast_category(typename std::iterator_traits<InputIt>::iterator_category{}, first, last, root);
	}
	template<class InputIt, typename Size>
	auto broadcast_n(InputIt first, Size count, int root = 0){
		return broadcast_n_category(typename detail::iterator_category<InputIt>{}, first, count, root);
	}
private:
	template<class RandomAccessIt>
	auto broadcast_category(std::random_access_iterator_tag, RandomAccessIt first, RandomAccessIt last, int root){
		broadcast_n(first, std::distance(first, last), root);
	}
	template<class ContiguousIt, typename Size>
	auto broadcast_n_category(detail::contiguous_iterator_tag, ContiguousIt first, Size count, int root){
		return broadcast_n_contiguous_builtinQ(detail::is_basic<typename std::iterator_traits<ContiguousIt>::value_type>{}, first, count, root);
	}
	template<class ContiguousIt, typename Size, class sendtype = detail::basic_datatype<typename std::iterator_traits<ContiguousIt>::value_type>>
	auto broadcast_n_contiguous_builtinQ(std::true_type, ContiguousIt first, Size count, int root){
		int s = MPI_Bcast( detail::data(first), count, sendtype::value, root, impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot broadcast");
	}
	template<class ContiguousIt, typename Size>
	void broadcast_n_contiguous_builtinQ(std::false_type, ContiguousIt first, Size count, int root);
#endif
	template<class It, typename Size>
	auto ibroadcast_n(
		It first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size count,
		int root
	){
		request r;
		int s = MPI_Ibcast(
			detail::data(first), count, 
			detail::basic_datatype<typename std::iterator_traits<It>::value_type>{},
			root, impl_, &r.impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot ibroadcast"};
		return r;
	}
	template<class It, typename Size>
	auto broadcast_n(
		It first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size count,
		int root
	){
		auto e = static_cast<enum error>(
			MPI_Bcast(
				detail::data(first), count, 
				detail::basic_datatype<typename std::iterator_traits<It>::value_type>{},
			root, impl_)
		);
		if(e != boost::mpi3::error::success) throw std::system_error{e, "cannot broadcast"};
		return first + count;
	}
	template<class It>
	auto ibroadcast(
		It first, It last, 
			detail::random_access_iterator_tag,
		int root
	){
		return ibroadcast_n(
			first, std::distance(first, last), root
		);
	}
	template<class It>
	auto broadcast(
		It first, It last, 
			detail::random_access_iterator_tag,
		int root
	){
		return broadcast_n(
			first, std::distance(first, last), root
		);
	}
public:
	template<class It, typename Size>
	auto ibroadcast_n(It first, Size count, int root = 0){
		return ibroadcast_n(
			first,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			count, 
			root
		);
	}
	template<class It, typename Size>
	auto broadcast_n(It first, Size count, int root = 0){
		return broadcast_n(
			first,
				detail::iterator_category_t<It>{},
				detail::value_category_t<typename std::iterator_traits<It>::value_type>{},
			count, 
			root
		);
	}
	template<class It>
	auto ibroadcast(It first, It last, int root = 0){
		return ibroadcast(
			first, last,
				detail::iterator_category_t<It>{},
			root
		);
	}
	template<class It>
	auto broadcast(It first, It last, int root = 0){
		return broadcast(
			first, last,
				detail::iterator_category_t<It>{},
			root
		);
	}
	template<class T, class Op = std::plus<> >
	void reduce_value(T const& t, T& ret, Op op = {}, int root = 0){
		reduce_n(std::addressof(t), 1, std::addressof(ret), op, root); 
	}
	template<class T, class Op = std::plus<> >
	auto reduce_value(T const& t, Op op = {}, int root = 0){
		T ret = T(0);
		reduce_value(t, ret, op, root); // if(rank() == root) return optional<T>(ret);
		return ret;
	}
//	template<class It1, class It2, class Op>
//	auto reduce(It1 first, It1 last, It2 d_first, Op op, int root){
//		return reduce_n(first, std::distance(first, last), d_first, op, root);
	//	return reduce_category(typename std::iterator_traits<It1>::iterator_category{}, first, last, d_first, op, root);
//	}
#if 0
	template<
		class It1, class Size, class It2, class Op, 
		class V1 = typename std::iterator_traits<It1>::value_type, class V2 = typename std::iterator_traits<It2>::value_type,
		class P1 = decltype(detail::data(It1{})), class P2 = decltype(detail::data(It2{})),
		typename = std::enable_if_t<std::is_same<V1, V2>{}>,
		class PredefinedOp = predefined_operation<Op>
	>
	auto reduce_n(It1 first, Size count, It2 d_first, Op /*op*/, int root){
		int s = MPI_Reduce(detail::data(first), detail::data(d_first), count, detail::basic_datatype<V1>{}, PredefinedOp{}, root, impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot reduce n");
	}
#endif
/*	template<class It1, class Size, class It2, class Op>
	auto reduce_n(It1 first, Size count, It2 d_first, Op op, int root = 0){
		return reduce_n(
			first, 
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			count, 
			d_first,
				detail::iterator_category_t<It2>{},
				detail::value_category_t<typename std::iterator_traits<It2>::value_type>{},
			op, root
		);			
	}*/
	template<class It1, class Size, class It2, class Op, class PredefinedOp>
	auto reduce_n(
		It1 first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size count, 
		It2 d_first,
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Op,
		PredefinedOp,
		int root
	){
		static_assert(std::is_same<typename std::iterator_traits<It1>::value_type, typename std::iterator_traits<It2>::value_type>{}, "!");
		int s = MPI_Reduce(
			detail::data(first), detail::data(d_first), count, 
			detail::basic_datatype<typename std::iterator_traits<It2>::value_type>{}, PredefinedOp{}, 
			root, impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot reduce n");
		return rank()==root?d_first + count:d_first;
	}
			
	template<class It1, class Size, class It2, class Op>
	auto reduce_n(It1 first, Size count, It2 d_first, Op op, int root = 0){
		return reduce_n(
			first, 
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			count, 
			d_first,
				detail::iterator_category_t<It2>{},
				detail::value_category_t<typename std::iterator_traits<It2>::value_type>{},
			op, 
				predefined_operation<Op>{},
			root
		);			
	}
	template<class T, class Op = std::plus<> >
	void all_reduce_value(T const& t, T& ret, Op op = {}){
		all_reduce_n(std::addressof(t), 1, std::addressof(ret), op); 
	}
	template<class T, class Op = std::plus<> >
	auto all_reduce_value(T const& t, Op op = {}){
		T ret = T(0);
		all_reduce_value(t, ret, op); // if(rank() == root) return optional<T>(ret);
		return ret;
	}
public:
	template<
		class It1, class Size, class It2, class Op = std::plus<>, 
		class V1 = typename std::iterator_traits<It1>::value_type, class V2 = typename std::iterator_traits<It2>::value_type,
		class P1 = decltype(detail::data(It1{})), class P2 = decltype(detail::data(It2{})),
		typename = std::enable_if_t<std::is_same<V1, V2>{}>,
		class PredefinedOp = predefined_operation<Op>,
		typename = decltype(predefined_operation<Op>::value)
	>
	auto all_reduce_n(It1 first, Size count, It2 d_first, Op /*op*/ = {}){
		using detail::data;
		int s = MPI_Allreduce(data(first), detail::data(d_first), count, detail::basic_datatype<V1>{}, PredefinedOp{}/*op*/, impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot reduce n");
	}
	template<class It1, class Size, class It2, class Op = std::plus<>, class PredefinedOp>
	auto all_reduce_n(
		It1 first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size count, 
		It2 d_first,
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Op,
			PredefinedOp
	){
		static_assert(std::is_same<typename std::iterator_traits<It1>::value_type, typename std::iterator_traits<It2>::value_type>{}, "!");
		using detail::data;
		int s = MPI_Allreduce(detail::data(first), detail::data(d_first), count, detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{}, PredefinedOp{}/*op*/, impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot reduce n"};
		return d_first + count;
	}
	template<class It1, class Size, class It2, class Op = std::plus<>>
	auto all_reduce_n(It1 first, Size count, It2 d_first, Op op = {}){
		return all_reduce_n(
			first, 
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			count, 
			d_first,
				detail::iterator_category_t<It2>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			op, 
				predefined_operation<Op>{}
		);
	}
	template<typename It1, typename It2, class Op = std::plus<>>
	auto all_reduce(It1 first, It1 last, It2 d_first, Op op = {}){
		return all_reduce_n(first, std::distance(first, last), d_first, op);
	}
	template<class T> static auto data_(T&& t){
		using detail::data;
		return data(std::forward<T>(t));
	}

	template<
		class It1, class Size, class Op, 
		class V1 = typename std::iterator_traits<It1>::value_type, class P1 = decltype(data_(It1{})), 
		class PredefinedOp = predefined_operation<Op>
	>
	auto all_reduce_in_place_n(It1 first, Size count, Op /*op*/){ // TODO check why op is not used 
		int s = MPI_Allreduce(MPI_IN_PLACE, data_(first), count, detail::basic_datatype<V1>{}, PredefinedOp{}, impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot all reduce n");
	}
	template<
		class It1, class Size, class Op, 
		class V1 = typename std::iterator_traits<It1>::value_type, class P1 = decltype(detail::data(It1{})), 
		class PredefinedOp = predefined_operation<Op>,
		typename = decltype(predefined_operation<Op>::value)
	>
	auto all_reduce_n(It1 first, Size count, Op op){
		return all_reduce_in_place_n(first, count, op);
	}

	template<
		class It1, class Size, class Op, 
		class V1 = typename std::iterator_traits<It1>::value_type, class P1 = decltype(data_(It1{})), 
		class PredefinedOp = predefined_operation<Op>
	>
	auto reduce_in_place_n(It1 first, Size count, Op /*op*/, int root = 0){
		int s = MPI_Reduce(MPI_IN_PLACE, data_(first), count, detail::basic_datatype<V1>{}, PredefinedOp{}, root, impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot reduce n"};
	}
	template<
		class It1, class Size, class Op = std::plus<>, 
		class V1 = typename std::iterator_traits<It1>::value_type, class P1 = decltype(detail::data(It1{})), 
		class PredefinedOp = predefined_operation<Op>
	>
	auto reduce_n(It1 first, Size count, Op op = {}, int root = 0){
		if(rank() == root) return reduce_in_place_n(first, count, op, root);
		return reduce_n(first, 0, nullptr, op, root);
	}
	template<
		class It1, class Op, 
		class V1 = typename std::iterator_traits<It1>::value_type, class P1 = decltype(detail::data(It1{})), 
		class PredefinedOp = predefined_operation<Op>
	>
	auto reduce_in_place(It1 first, It1 last, Op op, int root = 0){
		return reduce_in_place_n(first, std::distance(first, last), op, root);
	}
	template<
		class It1, class Op = std::plus<>, 
		class V1 = typename std::iterator_traits<It1>::value_type, class P1 = decltype(detail::data(It1{})), 
		class PredefinedOp = predefined_operation<Op>
	>
	auto reduce(It1 first, It1 last, Op op = {}, int root = 0){
		return reduce_n(first, std::distance(first, last), op, root);
	}
	
private:

	template<class ReducePolicy, class It1, class Size, class It2, class Op, 
		class V1 = typename std::iterator_traits<It1>::value_type, 
		class V2 = typename std::iterator_traits<It2>::value_type 
	>
	auto reduce_n(ReducePolicy rp, It1 first, Size count, It2 d_first, Op op, int root = 0){
		return reduce_n_dispatch(rp, 
			std::integral_constant<bool, detail::is_basic<V1>{} and detail::is_basic<V2>{} and detail::is_contiguous<It1>{} and detail::is_contiguous<It2>{}>{}, 
			first, count, d_first, op, root
		);
	}
	template<class ReducePolicy, class RandomAccessIt1, class It2, class Op>
	auto reduce_category(ReducePolicy rp, std::random_access_iterator_tag, RandomAccessIt1 first, RandomAccessIt1 last, It2 d_first, Op op, int root){
		return reduce_n(rp, first, std::distance(first, last), d_first, op, root);
	}
//	template<class It1, class Size, class It2>
//	auto reduce_n_category(detail::contiguous_iterator_tag, detail::contiguous_iterator_tag, It1 first, Size count, It2 d_first, Op op, int root){
//	}
	template<class It1, class Size, class It2, class Op, 
		class V1 = typename std::iterator_traits<It1>::value_type, 
		class V2 = typename std::iterator_traits<It2>::value_type 
	>
	auto reduce_n_dispatch(all_reduce_mode rp, std::true_type, It1 first, Size count, It2 d_first, Op op, int /*root*/ = 0){
		int s = rp(
			detail::data(first)  , 
			detail::data(d_first), count, detail::basic_datatype<V1>{},
			op.impl_, impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot scatter");
	}
	template<class It1, class Size, class It2, class Op,
		class V1 = typename std::iterator_traits<It1>::value_type, 
		class V2 = typename std::iterator_traits<It2>::value_type 
	>
	auto reduce_n_dispatch(reduce_mode rp, std::true_type, It1 first, Size count, It2 d_first, Op op, int root = 0){
		int s = rp(
			detail::data(first)  , 
			detail::data(d_first), count, detail::basic_datatype<V1>{},
			op.impl_,
			root, impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot scatter");
	}
	template<class It1, class Size, class It2, class Op,
		class V1 = typename std::iterator_traits<It1>::value_type, 
		class V2 = typename std::iterator_traits<It2>::value_type 
	>
	auto reduce_n_dispatch(ireduce_mode rp, std::true_type, It1 first, Size count, It2 d_first, Op op, int root = 0){
		request ret;
		int s = rp(
			detail::data(first)  , 
			detail::data(d_first), count, detail::basic_datatype<V1>{},
			op.impl_,
			root, impl_, &ret.impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot scatter");
		return ret;
	}

public:
	template<class InputIt1, class Size, class InputIt2, 
		class V1 = typename std::iterator_traits<InputIt1>::value_type, 
		class V2 = typename std::iterator_traits<InputIt2>::value_type 
	>
	void scatter_n(InputIt1 first, Size count, InputIt2 d_first, int root = 0){
		return scatter_n_dispatch( 
			std::integral_constant<bool, 
				detail::is_basic<V1>{} and detail::is_basic<V2>{}
				and detail::is_contiguous<InputIt1>::value and detail::is_contiguous<InputIt2>::value
			>{}, first, count, d_first, root
		);
	}
	template<class InputIt1, class Size, class InputIt2, 
		class V1 = typename std::iterator_traits<InputIt1>::value_type, 
		class V2 = typename std::iterator_traits<InputIt2>::value_type 
	>
	void scatter_from_n(InputIt2 d_first, Size count, InputIt1 first, int root = 0){
		return scatter_from_n_dispatch( 
			std::integral_constant<bool, 
				detail::is_basic<V1>{} and detail::is_basic<V2>{}
				and detail::is_contiguous<InputIt1>::value and detail::is_contiguous<InputIt2>::value
			>{}, d_first, count, first, root
		);
	}
	template<class InputIt1, class InputIt2>
	void scatter(InputIt1 first, InputIt1 last, InputIt2 d_first, int root = 0){
		return scatter_category(typename std::iterator_traits<InputIt1>::iterator_category{}, first, last, d_first, root);
	}
	template<class InputIt1, class InputIt2>
	void scatter_from(InputIt1 d_first, InputIt1 d_last, InputIt2 first, int root = 0){
		return scatter_from_category(typename std::iterator_traits<InputIt1>::iterator_category{}, d_first, d_last, first, root);
	}
	template<class InputIt1, class T>
	void scatter_value(InputIt1 first, T& t, int root = 0){
		return scatter_n(first, size(), std::addressof(t), root);
	}
	template<class InputIt1>
	typename std::iterator_traits<InputIt1>::value_type 
	scatter_value(InputIt1 first, int root = 0){
		typename std::iterator_traits<InputIt1>::value_type ret;
		scatter_value(first, ret, root);
		return ret;
	}
	template<class T>
	T scatter_value(std::vector<T> const& v, int root = 0){
		if(rank() == root) assert( static_cast<std::ptrdiff_t>(v.size()) == size() );
		return scatter_value(v.begin(), root);
	}
private:
	template<class InputIt1, class InputIt2>
	void scatter_category(std::random_access_iterator_tag, InputIt1 first, InputIt1 last, InputIt2 d_first, int root){
		return scatter_n(first, std::distance(first, last), d_first, root);
	}
	template<class InputIt1, class InputIt2>
	void scatter_from_category(std::random_access_iterator_tag, InputIt1 d_first, InputIt1 d_last, InputIt2 first, int root){
		return scatter_from_n(d_first, std::distance(d_first, d_last), first, root);
	}
	template<class ContIt1, class Size, class ContIt2,
		class V1 = typename std::iterator_traits<ContIt1>::value_type, 
		class V2 = typename std::iterator_traits<ContIt2>::value_type 
	>
	void scatter_n_dispatch(std::true_type, ContIt1 first, Size count, ContIt2 d_first, int root){
		if(count % size() != 0) throw std::runtime_error("not matching size for scatter, elements count is " + std::to_string(count) + ", comm size is " + std::to_string(size()));
		int s = MPI_Scatter(
			detail::data(first), count/size(), detail::basic_datatype<V1>{},
			detail::data(d_first), count/size(), detail::basic_datatype<V2>{},
			root, impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot scatter");
	}
	template<class ContIt1, class Size, class ContIt2,
		class V1 = typename std::iterator_traits<ContIt1>::value_type, 
		class V2 = typename std::iterator_traits<ContIt2>::value_type 
	>
	void scatter_from_n_dispatch(std::true_type, ContIt2 d_first, Size count, ContIt1 first,  int root){
		if(count % size() != 0) throw std::runtime_error("not matching size for scatter");
		int s = MPI_Scatter(
			detail::data(first), count, detail::basic_datatype<V1>{},
			detail::data(d_first), count, detail::basic_datatype<V2>{},
			root, impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot scatter");
	}

#if 0
	template<class IContiguousIterator, class OContiguousIterator, 
		class IT = typename std::iterator_traits<IContiguousIterator>::value_type,
		class OT = typename std::iterator_traits<OContiguousIterator>::value_type
	>
	void scatter_n(
		IContiguousIterator first , std::ptrdiff_t send_count, 
		OContiguousIterator result, std::ptrdiff_t receive_count, 
		int root = 0
	) const{
		// assumes: first, first + send_count is not aliased with result
		assert(send_count%size() == 0);
		assert(receive_count%size() == 0);
		MPI_Scatter(
			(void*)(&*first) , send_count/size()   , detail::basic_datatype<IT>::value, 
			(void*)(&*result), receive_count/size(), detail::basic_datatype<OT>::value, 
			root, impl_
		);
	}
#endif
//	template<class IContiguousIterator, class OContiguousIterator>
//	void scatter_n(IContiguousIterator first, std::ptrdiff_t send_count, OContiguousIterator result, int root = 0) const{
//		scatter_n(first, send_count, result, send_count, root);
//	}
public:
	template<class T, class It> 
	void all_gather_value(T const& t, It first){all_gather_n(std::addressof(t), 1, first);}
	template<class T> std::vector<T> 
	all_gather_value(T const& t){
		std::vector<T> ret(size());
		all_gather_value(t, ret.begin());
		return ret;
	}

//	template<class T> 
//	std::vector<T> gather_value(T const& t, int root = 0){
//		std::vector<T> ret((rank() == root)?size():0);
//		gather_value(t, ret.begin(), root);
//		return ret;
//	}
//	template<typename T, typename It> 
//	void gather_value(T const& t, It first, int root){
//		gather_n(std::addressof(t), 1, first, root);
//	}
	protected:
	template<class It, typename Size> 
	void advance(It& it, Size count){std::advance(it, count);}
	template<class It, typename Size> 
	void advance(It& it, Size s, int r){std::advance(it, rank()==r?s:0);}
	template<
		class GatherMode,
		typename It1, typename Size, typename It2, 
		typename = std::enable_if_t<
			std::is_same<
				typename std::iterator_traits<It1>::value_type, 
				typename std::iterator_traits<It2>::value_type
			>{}
		>, class... Root
	>
	auto a_gather_n(
		GatherMode gm,
		It1 first, 
			detail::contiguous_iterator_tag, 
			detail::memcopyable_tag,
		Size count, 
		It2 d_first, 
			detail::contiguous_iterator_tag,
			detail::memcopyable_tag,
		Root... root
	){
		int s = gm(
			detail::data(first), count*sizeof(*first), MPI_BYTE, 
			detail::data(d_first), count*sizeof(*d_first), MPI_BYTE, 
			root..., impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot gather");
		advance(d_first, count*size(), root...);
	//	std::advance(d_first, count);
		return d_first;
	}
	public:
/*	template<typename It1, typename Size, typename It2>
	auto gather_n(It1 first, Size count, It2 d_first, Size d_count, int root){
		return a_gather_n(
			gather_mode{},
			first, 
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			count, 
			d_first, 
				detail::iterator_category_t<It2>{},
				detail::value_category_t<typename std::iterator_traits<It2>::value_type>{},
			d_count,
			root
		);
	}*/
	template<typename It1, typename Size, typename It2>
	auto all_gather_n(
		It1 first  , Size count  , 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		It2 d_first, Size d_count,
			detail::contiguous_iterator_tag,
			detail::basic_tag
	){
		int s = MPI_Allgather(
			detail::data(first)  , count  , detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{},
			detail::data(d_first), d_count, detail::basic_datatype<typename std::iterator_traits<It2>::value_type>{},
			impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot scatter");
		return d_first + d_count*size();
	}
	template<typename It1, typename Size, typename It2>
	auto all_gather_n(It1 first, Size count, It2 d_first, Size d_count){
		return all_gather_n(
			first, count,
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			d_first, d_count,
				detail::iterator_category_t<It2>{},
				detail::value_category_t<typename std::iterator_traits<It2>::value_type>{}
		);
	}
	template<typename It1, typename Size, typename It2, typename CountsIt, typename DisplsIt>
	auto all_gatherv_n(It1 first, Size count, It2 d_first, CountsIt counts, DisplsIt displs){
		return all_gatherv_n(
			first, count, 
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			d_first,
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			counts, displs
		);
	}
	template<typename It1, typename Size, typename It2, typename CountsIt, typename DisplsIt>
	auto all_gatherv_n(
		It1 first, Size count, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		It2 d_first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		CountsIt counts, DisplsIt displs
	){
		int s = MPI_Allgatherv(
			detail::data(first), count,
			detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{},
			detail::data(d_first), detail::data(counts), detail::data(displs),
			detail::basic_datatype<typename std::iterator_traits<It2>::value_type>{},
			impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot all_gatherv");
		return d_first + detail::data(displs)[size()-1] + detail::data(counts)[size()-1];
	}
	template<typename It1, typename Size, typename It2>
	auto all_gather_n(
		It1 first  , Size count  , 
			detail::forward_iterator_tag,
			detail::value_unspecified_tag,
		It2 d_first, Size d_count,
			detail::forward_iterator_tag,
			detail::value_unspecified_tag
	){
		detail::package po(*this);
		package_oarchive poa(po);
		while(count--) poa << *(first++);
		int posize = po.size();
		std::vector<int> sizes(size());
		std::vector<int> displs(size());
		all_gather_n(&posize, 1, sizes.begin(), 1);
		partial_sum(sizes.begin(), sizes.end(), displs.begin()+1);
		detail::package pi(*this);
		int total = std::accumulate(sizes.begin(), sizes.end(), 0);
		pi.resize(total);
		all_gatherv_n(po.data(), po.size(), pi.data(), sizes.data(), displs.data());
		package_iarchive pia(pi);
		d_count*=size();
		while(d_count--) pia >> *(d_first++);
		return d_first;
	}
	template<typename It1, typename Size, typename It2>
	auto all_gather_n(It1 first, Size count, It2 d_first){
		return all_gather_n(first, count, d_first, count);
	}
	template<typename V, typename It2>
	auto all_gather_v(V const& v, It2 d_first){
		return all_gather_n(std::addressof(v), 1, d_first);
	}
	template<class Vector, typename V>
	auto all_gather_as(V const& v){
		Vector ret(size());
		all_gather_v(v, ret.begin());
		return ret;
	}
	template<typename It1, typename Size, typename It2>
	auto all_gatherv_n(It1 first, Size count, It2 d_first){
		std::vector<int> counts(size());
		std::vector<int> displs(size());
		int c = count;
		all_gather_n(&c, 1, counts.begin());
		partial_sum(counts.begin(), counts.end(), displs.begin()+1);
		return all_gatherv_n(first, count, d_first, counts.begin(), displs.begin());
	}
	template<class It1, typename Size1, class It2, class Size2>
	auto all_gather_n(
		It1 first,
			detail::forward_iterator_tag,
			detail::value_unspecified_tag,
		Size1 count, 
		It2 d_first,
			detail::forward_iterator_tag,
			detail::value_unspecified_tag,
		Size2 d_count
	){
		detail::package po(*this);
		package_oarchive poa(po);
		while(count--) poa << *(first++);
		int posize = po.size();
		std::vector<int> sizes(size());
		std::vector<int> displs(size());
		all_gather_n(&posize, 1, sizes.begin(), 1);
		partial_sum(sizes.begin(), sizes.end(), displs.begin()+1);
		detail::package pi(*this);
		int total = std::accumulate(sizes.begin(), sizes.end(), 0);
		pi.resize(total);
		all_gatherv_n(po.data(), po.size(), pi.data(), sizes.data(), displs.data());
		package_iarchive pia(pi);
		d_count*=size();
		while(d_count--) pia >> *(d_first++);
		return d_first;
	}

/*	template<typename It1, typename Size, typename It2>
	auto all_gather_n(It1 first, Size count, It2 d_first){
		return
			a_gather_n(
				all_gather_mode{},
				first, 
					detail::iterator_category_t<It1>{},
					detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
				count, 
				d_first, 
					detail::iterator_category_t<It2>{},
					detail::value_category_t<typename std::iterator_traits<It2>::value_type>{}
			);
	}*/
	public:
/*
	template<typename It1, typename It2>
	auto all_gather(It1 first, It1 last, It2 d_first){
		return
			a_gather(
				all_gather_mode{},
				first, last,
					detail::iterator_category_t<It1>{},
				d_first 
			);
	}*/
	protected:
/*
	template<typename GatherMode, typename It1, typename It2, typename... Root>
	auto a_gather(
		GatherMode gm, 
		It1 first, It1 last, 
			detail::random_access_iterator_tag,
		It2 d_first, Root... root
	){
		return a_gather_n(
			gm, 
			first, 
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			std::distance(first, last), 
			d_first, 
				detail::iterator_category_t<It2>{},
				detail::value_category_t<typename std::iterator_traits<It2>::value_type>{},
			root...
		);
	}
*/
//	auto gather_n_randomaccess_contiguous_builtin(std::true_type, std::true_type, ContiguousIt1 first, Size count, ContiguousIt2 d_first, Size d_count, int root){
//		int s = MPI_Gather(detail::data(first), count, datatype1{}, detail::data(d_first), d_count, datatype2{}, root, impl_);
//		if(s != MPI_SUCCESS) throw std::runtime_error("cannot gather");
//	}
public:
	template<class It1, typename Size1, class It2, class Itc, class Itd>
	auto gatherv_n(
		It1 first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size1 count,
		It2 d_first,
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Itc counts,
			detail::contiguous_iterator_tag,
		Itd displs,
			detail::contiguous_iterator_tag,
		int root
	){
		MPI_Gatherv(
			detail::data(first), count, detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{},
			detail::data(d_first), detail::data(counts), detail::data(displs), detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{},
			root, impl_
		);
	}
	template<class It1, typename Size1, class It2, class Itc, class Itd>
	auto gatherv_n(
		It1 first, Size1 count,
		It2 d_first, Itc counts, Itd displs,
		int root
	){
		return gatherv_n(
			first, 
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			count,
			d_first,
				detail::iterator_category_t<It2>{},
				detail::value_category_t<typename std::iterator_traits<It2>::value_type>{},
			counts, 
				detail::iterator_category_t<Itc>{},
			displs, 
				detail::iterator_category_t<Itd>{},
			root
		);
	}
	template<class It1, typename Size1, class It2, typename Size2>
	auto gather_n(
		It1 first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size1 count,
		It2 d_first,
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size2 d_count, 
		int root
	){
		int s = MPI_Gather(
			detail::data(first), count, 
			detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{},
			detail::data(d_first), d_count,
			detail::basic_datatype<typename std::iterator_traits<It2>::value_type>{},
			root, impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot gather");
		if(rank() == root) d_first += d_count*size();
		return d_first;
	}
	template<class It1, typename Size1, class It2, typename Size2>
	auto igather_n(
		It1 first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size1 count,
		It2 d_first,
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size2 d_count, 
		int root
	){
		request ret;
		int s = MPI_Igather(
			detail::data(first)  , count  , detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{},
			detail::data(d_first), d_count, detail::basic_datatype<typename std::iterator_traits<It2>::value_type>{}, 
			root, impl_, &ret.impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot Igather"};
		return ret;
	}
	template<class It1, typename Size1, class It2, typename Size2>
	auto iall_gather_n(
		It1 first, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size1 count,
		It2 d_first,
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		Size2 d_count
	){
		request ret;
		int s = MPI_Iallgather(
			detail::data(first)  , count  , detail::basic_datatype<typename std::iterator_traits<It1>::value_type>{},
			detail::data(d_first), d_count, detail::basic_datatype<typename std::iterator_traits<It2>::value_type>{}, 
			impl_, &ret.impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot Iallgather");
		return ret;
	}
	template<class It1, class Size1, class It2, class Size2>
	auto gather_n(It1 first, Size1 count, It2 d_first, Size2 d_count, int root){
		return gather_n(
			first, 
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			count, 
			d_first,
				detail::iterator_category_t<It2>{},
				detail::value_category_t<typename std::iterator_traits<It2>::value_type>{},
			d_count,
			root
		);
	}
	template<class It1, typename Size1, class It2, class Size2>
	auto gather_n(
		It1 first,
			detail::input_iterator_tag,
			detail::value_unspecified_tag,
		Size1 count, 
		It2 d_first,
			detail::input_iterator_tag,
			detail::value_unspecified_tag,
		Size2 d_count,
		int root
	){
		detail::package po(*this);
		package_oarchive poa(po);
		while(count--) poa << *(first++);
		int posize = po.size();
		std::vector<int> sizes(rank()==root?size():0);
		std::vector<int> displs(rank()==root?size():0);
		gather_n(&posize, 1, sizes.begin(), 1, root);
		partial_sum(sizes.begin(), sizes.end(), displs.begin()+1);
		detail::package pi(*this);
		int total = std::accumulate(sizes.begin(), sizes.end(), 0);
		pi.resize(total);
		gatherv_n(po.data(), po.size(), pi.data(), sizes.data(), displs.data(), root);
		if(rank() == root){
			package_iarchive pia(pi);
			d_count*=size();
			while(d_count--) pia >> *(d_first++);
		}
		return d_first;
	}
	template<class It1, class Size1, class It2, class Size2>
	auto igather_n(It1 first, Size1 count, It2 d_first, Size2 d_count, int root){
		return igather_n(
			first, 
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			count, 
			d_first,
				detail::iterator_category_t<It2>{},
				detail::value_category_t<typename std::iterator_traits<It2>::value_type>{},
			d_count,
			root
		);
	}
	template<class It1, class Size1, class It2, class Size2>
	auto iall_gather_n(It1 first, Size1 count, It2 d_first, Size2 d_count){
		return iall_gather_n(
			first, 
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			count, 
			d_first,
				detail::iterator_category_t<It2>{},
				detail::value_category_t<typename std::iterator_traits<It2>::value_type>{},
			d_count
		);
	}
	template<class It1, class Size1, class It2>
	auto gather_n(It1 first, Size1 count, It2 d_first, int root){
		return gather_n(first, count, d_first, count, root);
	}
	template<class It1, class Size1, class It2>
	auto igather_n(It1 first, Size1 count, It2 d_first, int root = 0){
		return igather_n(first, count, d_first, count, root);
	}
	template<class It1, class Size1, class It2>
	auto iall_gather_n(It1 first, Size1 count, It2 d_first){
		return iall_gather_n(first, count, d_first, count);
	}
	template<typename It1, typename It2>
	auto gather(
		It1 first, It1 last, 
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		It2 d_first, int root
	){
		return gather_n(first, std::distance(first, last), d_first, root);
	}
	template<typename It1, typename It2>
	auto gather(
		It1 first, It1 last, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		It2 d_first, int root
	){
		return gather_n(first, std::distance(first, last), d_first, root);
	}
	template<typename It1, typename It2>
	auto igather(
		It1 first, It1 last, 
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		It2 d_first, int root
	){
		return igather_n(first, std::distance(first, last), d_first, root);
	}
	template<typename It1, typename It2>
	auto iall_gather(
		It1 first, It1 last, 
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		It2 d_first
	){
		return iall_gather_n(first, std::distance(first, last), d_first);
	}
	template<class It1, class It2>
	auto gather(
		It1 first, It1 last, 
			detail::input_iterator_tag,
			detail::basic_tag,
		It2 d_first, int root
	){
		mpi3::vector<typename std::iterator_traits<It1>::value_type> buffer(first, last);
		return gather_n(buffer.data(), buffer.size(), d_first, root);
	}
	template<typename It1, typename It2>
	auto gather(It1 first, It1 last, It2 d_first, int root = 0){
		return gather(
			first, last,
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			d_first,
			root
		);
	}
	template<typename It1, typename It2>
	auto igather(It1 first, It1 last, It2 d_first, int root){
		return igather(
			first, last,
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			d_first,
			root
		);
	}
	template<typename It1, typename It2>
	auto iall_gather(It1 first, It1 last, It2 d_first){
		return iall_gather(
			first, last,
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			d_first
		);
	}
	template<class It1, class It2>
	auto gather(
		It1 first, It1 last, 
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		It2 d_first, It2 d_last,
			detail::contiguous_iterator_tag,
			detail::basic_tag,
		int root
	){
		return gather_n(
			first, std::distance(first, last), 
			d_first, std::distance(d_last, d_last), 
			root
		);
	}
	template<class It1, class It2>
	auto gather(
		It1 first, It1 last, 
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		It2 d_first, It2 d_last,
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		int root
	){
		return gather_n(
			first, std::distance(first, last), 
			d_first, std::distance(d_last, d_last), 
			root
		);
	}
	template<class It1, class It2>
	auto gather(
		It1 first, It1 last, 
			detail::forward_iterator_tag,
			detail::basic_tag,
		It2 d_first, It2 d_last,
			detail::random_access_iterator_tag,
			detail::value_unspecified_tag,
		int root
	){
		mpi3::vector<typename std::iterator_traits<It1>::value_type> v(first, last);
		return gather_n(
			v.data(), v.size(), 
			d_first, std::distance(d_first, d_last), 
			root
		);
	}
	public:
	template<typename It1, typename It2>
	auto gather(It1 first, It1 last, It2 d_first, It2 d_last, int root){
		return gather(
			first, last,
				detail::iterator_category_t<It1>{},
				detail::value_category_t<typename std::iterator_traits<It1>::value_type>{},
			d_first, d_last,
				detail::iterator_category_t<It2>{},
				detail::value_category_t<typename std::iterator_traits<It2>::value_type>{},
			root
		);
	}
//	template<class It1, class Size, class It2>
//	auto igather_n(It1 first, Size count, It2 d_first, int root = 0){return gather_n(igather_mode{}, first, count, d_first, root);}
/*
	template<class It1, class It2>
	//	[[nodiscard]] 
	auto igather(It1 first, It1 last, It2 d_first, int root = 0){
		return igather_category(
			typename std::iterator_traits<It1>::iterator_category{},
			first, last, d_first, root
		);
	}
	template<class It1, class It2>
	auto igather_category(std::random_access_iterator_tag, It1 first, It1 last, It2 d_first, int root = 0){
		return igather_n(first, std::distance(first, last), d_first, root);
	}*/
private:
/*	template<class GatherPolicy, class It1, class Size, class It2,
		class V1 = typename std::iterator_traits<It1>::value_type, 
		class V2 = typename std::iterator_traits<It2>::value_type 
	>
	auto gather_n(GatherPolicy gp, It1 first, Size count, It2 d_first, int root = 0){
		return gather_n_dispatch(gp, 
			std::integral_constant<bool, 
				detail::is_basic<V1>{} and detail::is_basic<V2>{}
				and detail::is_contiguous<It1>::value and detail::is_contiguous<It2>::value
			>{}, first, count, d_first, root);
	}*/
/*	template<class GatherPolicy, class It1, class It2>
	auto gather(GatherPolicy gp, It1 first, It1 last, It2 d_first, int root){
		return gather_category(gp, typename std::iterator_traits<It1>::iterator_category{}, first, last, d_first, root);
	}*/
/*	template<class GatherPolicy, class RandomAccessIt1, class It2>
	auto gather_category(GatherPolicy gp, std::random_access_iterator_tag, RandomAccessIt1 first, RandomAccessIt1 last, It2 d_first, int root){
		return gather_n(gp, first, std::distance(first, last), d_first, root);
	}*/
	template<class It1, class Size, class It2,
		class V1 = typename std::iterator_traits<It1>::value_type, 
		class V2 = typename std::iterator_traits<It2>::value_type 
	>
	auto gather_n_dispatch(gather_mode g, std::true_type, It1 first, Size count, It2 d_first, int root = 0){
		int s = g(
			detail::data(first)  , count, detail::basic_datatype<V1>{},
			detail::data(d_first), count, detail::basic_datatype<V2>{},
			root, impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot scatter");
	}
	template<class It1, class Size, class It2,
		class V1 = typename std::iterator_traits<It1>::value_type, 
		class V2 = typename std::iterator_traits<It2>::value_type 
	>
	auto gather_n_dispatch(gather_mode g, std::false_type, It1 first, Size count, It2 d_first, int root = 0){
		int s = g(
			detail::data(first)  , count, detail::basic_datatype<V1>{},
			detail::data(d_first), count, detail::basic_datatype<V2>{},
			root, impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot scatter");
	}
	template<class It1, class Size, class It2,
		class V1 = typename std::iterator_traits<It1>::value_type, 
		class V2 = typename std::iterator_traits<It2>::value_type 
	>
	auto gather_n_dispatch(igather_mode g, std::true_type, It1 first, Size count, It2 d_first, int root = 0){
		request r;
		int s = g(
			detail::data(first)  , count, detail::basic_datatype<V1>{},
			detail::data(d_first), count, detail::basic_datatype<V2>{},
			root, impl_, &r.impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot scatter");
		return r;
	}
	template<class It1, class Size, class It2,
		class V1 = typename std::iterator_traits<It1>::value_type, 
		class V2 = typename std::iterator_traits<It2>::value_type 
	>
	auto gather_n_dispatch(all_gather_mode g, std::true_type, It1 first, Size count, It2 d_first, int /*root*/= 0){
		int s = g(
			detail::data(first)  , count, detail::basic_datatype<V1>{},
			detail::data(d_first), count, detail::basic_datatype<V2>{},
			impl_
		);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot scatter");
	}

//	template<class IContiguousIterator, class OContiguousIterator>
//	void gather_n(IContiguousIterator first, std::ptrdiff_t send_count, OContiguousIterator result, int root = 0){
//		gather_n(first, send_count, result, send_count, root);
//	}
//	template<class IContiguousIterator, class OContiguousIterator>
//	void scatter(IContiguousIterator first, IContiguousIterator last, OContiguousIterator result, int root = 0){
//		scatter_n(first, std::distance(first, last), result, root);
//	}
//	template<class IContiguousIterator, class OContiguousIterator>
//	void gather (IContiguousIterator first, IContiguousIterator last, OContiguousIterator result, int root = 0){
//		gather_n(first, std::distance(first, last), result, root);
//	}
public:
	std::string get_name() const{
		int len;
		char comm_name[MPI_MAX_OBJECT_NAME];
		int status = MPI_Comm_get_name(impl_, comm_name, &len);
		if(status != MPI_SUCCESS) throw std::runtime_error("cannot set name");
		return std::string(comm_name, len);
	}
	void set_name(std::string const& s){
		int status = MPI_Comm_set_name(impl_, s.c_str());
		if(status != MPI_SUCCESS) throw std::runtime_error("cannot get name");
	}
	std::string name() const{
		return get_name();
	}
	void name(std::string const& s){set_name(s);}

	communicator parent() const{
		communicator ret;
		MPI_Comm_get_parent(&ret.impl_);
		return ret;
	}
	communicator spawn(std::string const& argv0, int np) const{
		communicator intercomm;
		MPI_Comm_spawn(argv0.data(), MPI_ARGV_NULL, np, MPI_INFO_NULL, 0, MPI_COMM_SELF, &intercomm.impl_, MPI_ERRCODES_IGNORE );
		return intercomm;
	}
	communicator intercommunicator_create(int local_leader, communicator const& peer, int remote_leader, int tag = 0) const{
		communicator ret;
		int s = MPI_Intercomm_create(impl_, local_leader, peer.impl_, remote_leader, tag, &ret.impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot create intercommunicator");
		return ret;
	}
	communicator create(int local_leader, communicator const& peer, int remote_leader, int tag = 0) const{
		return intercommunicator_create(local_leader, peer, remote_leader, tag);
	}
	communicator create(group const& g) const;
	communicator create_group(group const& g, int tag) const;
	FILE* fopen(const char* filename, int amode = MPI_MODE_RDWR | MPI_MODE_CREATE);

inline std::string const& name(communicator::topology const& t){
	static std::map<communicator::topology, std::string> const names = {
		{communicator::topology::undefined, "undefined"}, 
		{communicator::topology::graph, "graph"},
		{communicator::topology::cartesian, "cartesian"}};
	return names.find(t)->second;
}

template<class T>
friend T operator+=(communicator& comm, T const& t){
	T val = comm.all_reduce_value(t, std::plus<>{});
	return val;
}
template<class T>
friend communicator& operator<<(communicator& comm, T const& t){
	comm.send_value(t);
	return comm;
}


};

struct strided_range{
	int first;
	int last;
	int stride = 1;
	int front() const{return first;}
	int back() const{return last - 1;}
	int size() const{return (last - first) / stride;}
};

#if 0
struct group{
	MPI_Group impl_;
//	static group empty(){return group();}
	group() : impl_{MPI_GROUP_EMPTY}{}
	group(communicator const& comm){MPI_Comm_group(&comm, &impl_);}
//	template<class ContiguousIntIterator>
//	group(group const& other, ContiguousIntIterator ranks_begin, std::size_t n){
//		MPI_Group_incl(other.impl_, n, ranks_begin, &impl_);
//	}
	group(group const& other, std::vector<int> ranks){// : group(other, ranks.data(), ranks.size()){
		MPI_Group_incl(other.impl_, ranks.size(), ranks.data(), &impl_);
	}
	~group(){MPI_Group_free(&impl_);}
	int rank() const{
		int rank = -1; 
		int s = MPI_Group_rank(impl_, &rank);
		if(s != MPI_SUCCESS) throw std::runtime_error("rank not available");
		return rank;
	}
//	bool is_null() const{return MPI_GROUP_NULL == impl_;}
//	operator bool() const{return not is_null();}
	int size() const{
		int size = -1;
		MPI_Group_size(impl_, &size); 
		return size;
	}
	bool empty() const{
		return size() == 0;
	}

	mpi3::equality compare(group const& other) const{
		int result;
		MPI_Group_compare(impl_, other.impl_, &result);
		return static_cast<boost::mpi3::equality>(result);
	}
	friend mpi3::equality compare(group const& self, group const& other){return self.compare(other);}

	bool operator==(group const& other) const{return compare(other) == boost::mpi3::identical;}
	bool operator!=(group const& other) const{return not operator==(other);} 

	std::vector<int> translate_ranks(std::vector<int> const& ranks1, group const& other) const{
		std::vector<int> ranks2(ranks1.size());
		MPI_Group_translate_ranks(impl_, ranks1.size(), ranks1.data(), other.impl_, ranks2.data());
		return ranks2;
	}
	std::vector<int> translate_ranks(group const& other) const{
		std::vector<int> ranks1(other.size());
		for(int i = 0; i != (int)ranks1.size(); ++i) ranks1[i] = i;
	//	std::iota(ranks1.begin(), ranks1.end(), 0);
		return translate_ranks(ranks1, other);
	}
public:
	template<class Iterator, class Size>
	group include_n(Iterator it, Size count) const{
		return include_n_dispatch(typename std::iterator_traits<Iterator>::value_type{}, detail::iterator_category<Iterator>{}, it, count); 
	}
	template<class Iterator>
	group include(Iterator first, Iterator last) const{return include_n(first, std::distance(first, last));}
	template<class Range>
	group include(Range const& r) const{using std::begin; using std::end; return include(begin(r), end(r));}
private:
	template<class Iterator, class Size>
	group include_n_dispatch(int, detail::contiguous_iterator_tag, Iterator first, Size count) const{
		group ret;
		MPI_Group_incl(impl_, count, detail::data(first), &ret.impl_);
		return ret;
	}
	template<class Iterator, class Size, class Value, class = decltype(int(Value{}))>
	group include_n_dispatch(Value, std::input_iterator_tag, Iterator first, Size count) const{
		std::vector<int> v(count); std::copy_n(first, count, v.begin());
		return include_n(v.data(), v.size()); 
	}
public:
	group include(std::vector<int> const& ranks){return group(*this, ranks);}
	template<class ContiguousIntIterator, class Size>
	group exclude_n(ContiguousIntIterator it, Size count) const{
		group ret;
		MPI_Group_excl(impl_, count, boost::mpi3::detail::data(it), &ret.impl_);
		return ret;
	}
	template<class ContiguousIntIterator>
	group exclude(ContiguousIntIterator first, ContiguousIntIterator last) const{
		return exclude_n(first, std::distance(first, last));
	}
	template<class ContiguousRangeIterator, class Size>
	group range_exclude_n(ContiguousRangeIterator it, Size n) const{
		group ret;
		using int3 = int[3];
		auto ptr = const_cast<int3*>(reinterpret_cast<int3 const*>(boost::mpi3::detail::data(it)));
		int status = MPI_Group_range_excl(impl_, n, ptr, &ret.impl_);
		if(status != MPI_SUCCESS) throw std::runtime_error("cannot exclude");
		return ret;
	}
	template<class ContiguousRangeIterator, class Size>
	group range_include_n(ContiguousRangeIterator it, Size n) const{
		group ret;
		using int3 = int[3];
		auto ptr = const_cast<int3*>(reinterpret_cast<int3 const*>(boost::mpi3::detail::data(it)));
		int status = MPI_Group_range_incl(impl_, n, ptr, &ret.impl_);
		if(status  != MPI_SUCCESS) throw std::runtime_error("cannot include");
		return ret;
	}

	template<class ContiguousRangeIterator>
	group range_exclude(ContiguousRangeIterator first, ContiguousRangeIterator last) const{
		return range_exclude_n(first, std::distance(first, last));
	}
	group range_exclude(std::vector<boost::mpi3::strided_range> const& ranges){
		return range_exclude(ranges.begin(), ranges.end());
	}
	template<class ContiguousRangeIterator>
	group range_include(ContiguousRangeIterator first, ContiguousRangeIterator last) const{
		return range_include_n(first, std::distance(first, last));
	}
	group range_include(std::vector<boost::mpi3::strided_range> const& ranges){
		return range_include(ranges.begin(), ranges.end());
	}
	group set_union(group const& other) const{
		group ret;
		MPI_Group_union(impl_, other.impl_, &ret.impl_); 
		return ret;
	}
	friend group set_union(group const& self, group const& other){
		return self.set_union(other);
	}
	group set_intersection(group const& other) const{
		group ret;
		MPI_Group_intersection(impl_, other.impl_, &ret.impl_);
		return ret;
	}
	group intersection(group const& other) const{return set_intersection(other);}
	friend group set_intersection(group const& self, group const& other){
		return self.set_intersection(other);
	}
	friend group intersection(group const& self, group const& other){
		return self.intersection(other);
	}
//	template<class Iterator, class Size>
//	group include_n(Iterator it, Size n) const{
//		group ret;
//		MPI_Group_incl(impl_, n, std::addressof(*it), &ret.impl_);
//		return ret;
//	}
#if 0
	template<class ContiguousIntIterator>
	group include(ContiguousIntIterator it1, ContiguousIntIterator it2) const{
		return include_n(it1, std::distance(it1, it2));
	}
	template<class Range>
	group include(Range const& r) const{
		using std::begin;
		using std::end;
		return include(begin(r), end(r));
	}
#endif
};
#endif

//inline communicator::communicator(communicator const& other, struct group const& g, int tag) : communicator(){
//	MPI_Comm_create_group(other.impl_, g.impl_, tag, &impl_);
//}

//package communicator::make_package(int n){return package(*this, n);}

/*
inline void communicator::send_value(package const& p, int dest, int tag){
	int s = MPI_Send(p.data(), p.buffer_.size(), MPI_PACKED, dest, tag, impl_);
	if(s != MPI_SUCCESS) throw std::runtime_error("cannotsend_value package");
}*/


/*
template<class CommunicationMode, class ContiguousIterator, class Size, class ValueType = typename std::iterator_traits<ContiguousIterator>::value_type>
void communicator::send_n_randomaccess_contiguous_builtin(CommunicationMode cm, blocking_mode, std::false_type, ContiguousIterator first, Size n, int dest, int tag){
	package p(*this);
	detail::package_oarchive poa(p);
	while( n ){
		poa << *first;
		++first;
		--n;
	}
	p.send(dest, tag);
//	return first
}*/

/*
template<class CommunicationMode, class ContiguousIterator, class Size, class ValueType = typename std::iterator_traits<ContiguousIterator>::value_type>
void communicator::receive_n_randomaccess_contiguous_builtin(CommunicationMode cm, blocking_mode, std::false_type, ContiguousIterator first, Size n, int dest, int tag){
//	assert(0);
	package p(*this);
	p.receive(dest, tag);
	detail::package_iarchive pia(p);
	while(n){
		pia >> *first;
		++first;
		--n;
	}
}*/

/*
template<class ContiguousIt, typename Size>
void communicator::broadcast_n_contiguous_builtinQ(std::false_type, ContiguousIt first, Size count, int root){
	package p(*this);
	if(rank() == root){
		detail::package_oarchive poa(p);
		while(count){
			poa << *first;
			++first;
			--count;
		}
	}
	p.broadcast(root);
	if(rank() != root){
		detail::package_iarchive pia(p);
		while(count){
			pia >> *first;
			++first;
			--count;
		}
	}
}
*/

}}

//BOOST_SERIALIZATION_REGISTER_ARCHIVE(boost::mpi3::package_oarchive)
//BOOST_SERIALIZATION_USE_ARRAY_OPTIMIZATION(boost::mpi3::detail::package_oarchive)
//BOOST_SERIALIZATION_USE_ARRAY_OPTIMIZATION(boost::mpi3::detail::package_iarchive)

#ifdef _TEST_MPI3_COMMUNICATOR

#include "../mpi3/main.hpp"
#include "../mpi3/version.hpp"

#include<iostream>

using std::cout;
namespace mpi3 = boost::mpi3;

class V{
	mpi3::communicator comm_;
	public:
	V(mpi3::communicator const& c) : comm_(c){}
	V(mpi3::communicator&& c) : comm_(std::move(c)){}
};

int mpi3::main(int, char*[], mpi3::communicator world){

	static_assert(std::is_nothrow_constructible<mpi3::communicator>::value, "MyType should be noexcept MoveConstructible");

//	auto worldcopy1 = world;
//	auto worldcopy2 = std::move(worldcopy1);
//	V v(worldcopy);
//	V v2(std::move(v));

	if(world.rank() == 0) cout << "MPI version " <<  mpi3::version() << '\n';
//	if(world.rank() == 0) cout << "Topology: " << name(world.topo()) << '\n';

	cout << "MPI_ERR_COMM = " << MPI_ERR_COMM << '\n';

	mpi3::communicator comm;
	assert(!comm);
//	cout << comm.rank() << '\n';

	mpi3::communicator comm2 = world;
	assert(comm2);
	assert(comm2.size() == world.size());
	assert(comm2 == world);
	assert(&comm2 != &world);

	mpi3::communicator comm3 = world;//.duplicate();
	assert(comm3);
	assert(comm3 == world);
	assert(&comm3 != &world);
	comm = comm2;
	assert(&comm != &comm2);

//	world2 = world;

	return 0;
#if 0
//	boost::mpi3::communicator newcomm = world;
	{
		int color = world.rank()/3;
		communicator row_comm;
		row_comm = world.split(color);
		world.barrier();
		std::cout << std::to_string(world.rank()) + " " + std::to_string(row_comm.rank()) + "\n";// << std::endl;
		world.barrier();
	}
	{
		communicator row_comm = world/3;
		world.barrier();
		std::cout << std::to_string(world.rank()) + " " + std::to_string(row_comm.rank()) + "\n";// << std::endl;
		world.barrier();
	}

	world.barrier();
	if(world.rank() == 0) cout << "prime communicator" << '\n';
	world.barrier();

	{
	//	group world_group(world);
	//	const int ranks[4] = {2, 3, 5, 7};
	//	group prime = world_group.include(ranks, ranks + 4);
	//	communicator prime_comm(world, prime);
		auto prime_comm = world.subcomm({2,3,5,7});
		cout << world.rank() << " -> " << prime_comm.rank() << "/" << prime_comm.size() << '\n';
#if 0
		if(communicator::null != prime_comm){
			cout << world.rank() << " -> " << prime_comm.rank() << "/" << prime_comm.size() << '\n';
		}else{
			cout << world.rank() << " not in prime comm\n";
		}
#endif
	}

	world.barrier();
	if(world.rank() == 0) cout << "prime communicator" << '\n';
	world.barrier();

	if(0){
		auto prime = world.subcomm({2,3,5,7});
		if(prime.is_empty()){
	//	if (communicator::null != prime){
			cout << world.rank() << " -> " << prime.rank() << "/" << prime.size() << '\n';
		}else{
			cout << world.rank() << " not in prime comm\n";
		}
	}
#endif
}

#endif
#endif

