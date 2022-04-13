// UnorderedMap.cpp : Defines the entry point for the console application.
//
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <cmath>

  struct polymer_hash {
    std::size_t operator()(std::string const& s) const noexcept {
      std::size_t sum = 0;
      int r = s.size();
      int n = r - 1;
      int num = 0;
      for (int i = 0; i <= n; i++) {
        switch (std::tolower(s[i])) {
          case 'a':
            num = 0;
            break;
          case 'c':
            num = 1;
            break;
          case 't':
            num = 2;
            break;
          case 'g':
            num = 3;
            break;
        }
        sum += num * std::pow(r, n - i);
      }
      return sum;
    }
  };

class UnorderedMapPool;

bool operator==( const UnorderedMapPool& lhs,
                 const UnorderedMapPool& rhs );

bool operator!=( const UnorderedMapPool& lhs,
                 const UnorderedMapPool& rhs );

class UnorderedMapPool {
public:
  typedef std::string                       Key;
  typedef std::size_t                       Value;
  typedef Key                               key_type;
  typedef std::pair<Key const, Value>       value_type;
  typedef std::size_t                       size_type;
  typedef std::ptrdiff_t                    difference_type;
  typedef polymer_hash                      Hash;
  typedef Hash                              hash_function;
  typedef std::equal_to<Key>                KeyEqual;
  typedef value_type&                       reference;
  typedef const value_type&                 const_reference;
  class                                     iterator;
  class                                     const_iterator;

  struct node {
    node(value_type const& data = value_type(), node* next = NULL, node* iter_next = NULL)
      : data_(data),
        next_(next),
        iter_next_(iter_next)
    { }

    Key   const& key  () const { return data_.first;  }
    Value      & value()       { return data_.second; }
    Value const& value() const { return data_.second; }

    value_type data() { return data_; }
    node* next()      { return next_; }

    friend class UnorderedMapPool;
  private:
    value_type data_;
    node* next_;
    node* iter_next_;
  };

  UnorderedMapPool(size_type bucket_count = 50);
  UnorderedMapPool(const UnorderedMapPool& other); 

  UnorderedMapPool& operator=(const UnorderedMapPool& other ); 
  ~UnorderedMapPool() { clear(); delete[] buckets_; } 

  iterator begin(); 
  const_iterator begin() const; 
  const_iterator cbegin() const; 

  iterator end(); 
  const_iterator end() const; 
  const_iterator cend() const; 

  bool empty() const; 
  std::size_t size() const; 
  std::size_t max_size() const; 
  Value& at(const Key&);
  Value const& at(const Key&) const;

  void clear(); 
  std::pair<iterator, bool> insert(const value_type& value); 
  iterator insert(const_iterator hint, const value_type& value); 

  float max_load_factor() const; 
  void max_load_factor(float ml);
  float load_factor() const; 

  auto erase(const_iterator pos); 
  auto erase(const_iterator first, const_iterator last);  
  iterator InternalInsert(node** bucket, size_type pos, node* n);

  void swap(UnorderedMapPool& other); 
  Value& operator[](const Key& key); 
  std::size_t count(const Key& key) const; 
  iterator find(const Key& key);  
  const_iterator find(const Key& key) const; 
  iterator find_hint(size_t, Key const&);
  const_iterator find_hint(size_t, Key const&) const;
  std::pair<iterator, iterator> equal_range(const Key& key); 
  std::pair<const_iterator, const_iterator> equal_range(const Key& key) const; 

  void rehash(size_type count); 
  void reserve(size_type count); 
  auto key_eq() const; 

  friend bool operator==( const UnorderedMapPool& lhs,
                          const UnorderedMapPool& rhs ); 

  friend bool operator!=( const UnorderedMapPool& lhs,
                          const UnorderedMapPool& rhs ); 
  std::size_t bucket_count() const; 
  auto bucket_size(size_type n) const; 
public: // Set back to private
  node** buckets_;
  Hash hash_;
  KeyEqual pred_;
  std::size_t bucket_count_;
  size_type size_;
  node* begin_node_;
  float mlf_; // max load factor

  size_type bucket_hash(Key const& key) const {
    return hash_function()(key) % bucket_count();
  }
private:
  void init(size_type bucket_count) {
    buckets_ = new node*[bucket_count]();
    hash_ = polymer_hash();
    pred_ = KeyEqual();
    bucket_count_ = bucket_count;
    size_ = 0;
    begin_node_ = NULL;
    mlf_ = 0.0f;
  }

  float threshold() const { return bucket_count() * max_load_factor(); }

  node** Detach(node** start, node* node::*next, const key_type& key);

  void push_front(node*& head, node* newNode, node* node::*next) {
    newNode->*next = head;
    head = newNode;
  }

  void push_front_bucket_node(node*& head, node* newNode)
  { push_front(head, newNode, &node::next_); }
  void push_front_iterator_node(node*& head, node* newNode)
  { push_front(head, newNode, &node::iter_next_); }

  node*& bucket_at(Key const& key) const {
    return buckets_[bucket_hash(key)];
  }

 
};

UnorderedMapPool::UnorderedMapPool(size_type sz) {
  init(sz);
}

void UnorderedMapPool::reserve(std::size_t sz) {
  buckets_ = new node*[sz]();
}





auto UnorderedMapPool::bucket_size(size_type n) const {
  std::size_t count = 0;
  for (node* pos = buckets_[n]; pos != NULL; pos = pos->next()) {
    ++count;
  }
  return count;
}

bool UnorderedMapPool::empty() const {
  return size_ == 0;
}

float UnorderedMapPool::load_factor() const
{ return .75; }

std::size_t UnorderedMapPool::size() const
{ return size_; }

std::size_t UnorderedMapPool::bucket_count() const
{ return bucket_count_; }

std::size_t UnorderedMapPool::max_size() const
{ return std::numeric_limits<size_type>::max(); }

auto UnorderedMapPool::key_eq() const
{ return pred_; }

void UnorderedMapPool::clear() {
  for (std::size_t i = 0; i < bucket_count(); ++i) {
    for (node* p = buckets_[i]; p != NULL; ) {
      node* p_next = p->next();
      delete (p);
      p = p_next;
    }
    buckets_[i] = NULL;
  }
  begin_node_ = NULL;
  size_ = 0;
}

std::size_t UnorderedMapPool::count(const Key& key) const {
  node* pos = bucket_at(key);

  if (!pos) {
    return 0;
  }

  size_type c = 0;
  for (; pos; pos = pos->next()) {
    if (key_eq()(pos->key(), key)) {
      ++c;
    }
  }
  return c;
}


// Const-version
class UnorderedMapPool::const_iterator : public std::iterator<std::forward_iterator_tag, std::pair<Key const, Value>> {
public:
  explicit const_iterator(node* begin) : node_(begin) { }

  const_iterator& operator++() {
    node_ = node_->iter_next_;
    return *this;
  }

  const_iterator operator++(int) {
    const_iterator copy(*this);
    ++*this;
    return copy;
  }

  const_iterator operator+(size_type c) const {
    return const_iterator(*this) += c;
  }

  const_iterator& operator+=(size_type c) {
    while (c--) {
      ++*this;
    }
    return *this;
  }

  value_type const& operator*() const
  { return node_->data_; }

  value_type const* operator->() const
  { return &**this; }

  friend class UnorderedMapPool;
  friend class iterator;
  
  friend bool operator==(const const_iterator& lhs, const const_iterator& rhs)
  { return lhs.node_ == rhs.node_; }
  friend bool operator!=(const const_iterator& lhs, const const_iterator& rhs)
  { return !(lhs == rhs); }
private:
  node* node_;
};

// Non-const iterator
class UnorderedMapPool::iterator : public std::iterator<std::forward_iterator_tag, std::pair<Key const, Value>> {
public:
  explicit iterator(node* begin) : node_(begin) { }
  iterator() : node_(NULL) { }

  iterator& operator++() {
    node_ = node_->iter_next_;
    return *this;
  }

  iterator operator++(int) {
    iterator copy(*this);
    ++*this;
    return copy;
  }

  iterator operator+(size_type c) const {
    return iterator(*this) += c;
  }

  iterator& operator+=(size_type c) {
    while (c--) {
      ++*this;
    }
    return *this;
  }

  value_type const& operator*() const
  { return node_->data_; }

  value_type& operator*()
  { return node_->data_; }

  value_type const* operator->() const
  { return &(**this); }

  value_type* operator->()
  { return &(**this); }

  operator const_iterator() { return const_iterator(node_); }

  friend class UnorderedMapPool;
  friend class const_iterator;

  friend bool operator==(const const_iterator& lhs, const const_iterator& rhs);
  friend bool operator!=(const const_iterator& lhs, const const_iterator& rhs);
private:
  node* node_;
};

 auto UnorderedMapPool::InternalInsert(node** bucket, size_type pos, node* n) -> iterator {
    node*& head = bucket[pos];
    push_front_bucket_node(head, n);
    push_front_iterator_node(begin_node_, head);
    return iterator(head);
  }

  auto UnorderedMapPool::find_hint(size_type idx, Key const& key) -> iterator {
    node* dest = buckets_[idx];
    for (; dest != NULL; dest = dest->next_) {
      if (key_eq()(dest->key(), key)) {
        break;
      }
    }
    return iterator(dest);
  }
  auto UnorderedMapPool::find_hint(size_type idx, Key const& key) const -> const_iterator {
      return const_cast<UnorderedMapPool&>(*this).find_hint(idx, key);
  }

  auto UnorderedMapPool::Detach(node** start, node* node::*next, const key_type& key) -> node**
  {
    node** current = start;
    node* previous = NULL;
    for (; *current != NULL && !key_eq()(current[0]->key(), key); current = &(current[0]->*next)) {
      previous = *current;
    }

    if (*current) {
      if (previous) {
        previous->*next = current[0]->*next;
      }
    }
    return current;
  }

auto
UnorderedMapPool::insert(const_iterator /* hint */, const value_type& value) -> iterator {
  return insert(value).first;
}

auto UnorderedMapPool::find(Key const& key) -> iterator {
  return find_hint(bucket_hash(key), key);
}

auto UnorderedMapPool::find(Key const& key) const -> const_iterator
{ return const_cast<UnorderedMapPool&>(*this).find(key); }

auto UnorderedMapPool::begin() const -> const_iterator {
  // Find the first bucket that is non-null. Return an iterator to that node. It's that simple!
  return const_iterator(begin_node_);
}

auto UnorderedMapPool::begin() -> iterator {
  return iterator(begin_node_);
}

auto UnorderedMapPool::cbegin() const -> const_iterator
{ return begin(); }

auto UnorderedMapPool::end() -> iterator
{ return iterator(NULL); }

auto UnorderedMapPool::end() const -> const_iterator
{ return const_iterator(NULL); }

auto UnorderedMapPool::cend() const -> const_iterator
{ return end(); }

auto
UnorderedMapPool::equal_range(const Key& key) -> std::pair<iterator, iterator> {
  iterator it = find(key);
  return std::make_pair(it, it == end() ? it : it + 1);
}

auto
UnorderedMapPool::equal_range(const Key& key) const -> std::pair<const_iterator, const_iterator> {
  return const_cast<UnorderedMapPool&>(*this).equal_range(key);
}

void UnorderedMapPool::swap(UnorderedMapPool& other) {
  using std::swap;
  for (size_type i = 0; i < std::min(bucket_count(), other.bucket_count()); ++i) {
    swap(buckets_[i], other.buckets_[i]);
  }
  swap(buckets_, other.buckets_);
  swap(size_, other.size_);
  swap(mlf_, other.mlf_);
  swap(hash_, other.hash_);
  swap(begin_node_, other.begin_node_);
}

bool operator==( const UnorderedMapPool& lhs,
                 const UnorderedMapPool& rhs )
{
  if (lhs.size() != rhs.size()) return false;
  return std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

bool operator!=( const UnorderedMapPool& lhs,
                 const UnorderedMapPool& rhs )
{
  return !(lhs == rhs);
}

float UnorderedMapPool::max_load_factor() const {
  if (size() == 0) return 1.0f;
  return (float)bucket_count() / (float)size();
}

void UnorderedMapPool::max_load_factor(float ml) {
  mlf_ = ml;
}

auto UnorderedMapPool::at(const Key& key) -> Value& {
  size_type bucket_pos = bucket_hash(key);
  iterator it = find_hint(bucket_pos, key);
  if (it == end()) {
    throw std::invalid_argument("Could not find key");
  }
  return it->second;
}

auto UnorderedMapPool::at(const Key& key) const -> Value const& {
  return const_cast<UnorderedMapPool&>(*this).at(key);
}
#include<iostream>
using namespace std;
UnorderedMapPool& UnorderedMapPool::operator=(const UnorderedMapPool& other ) { // copy
  //rehash(other.bucket_count());

  delete[] buckets_;
  buckets_ = new node*[other.bucket_count()];
  //
  int c =0;
  for (const_iterator it = other.begin(); it != other.end(); ++it) {
    insert(*it);
  }
  bucket_count_ = other.bucket_count_;
  size_ = other.size_;
  begin_node_ = other.begin_node_;
  mlf_ = other.mlf_;
  hash_ = other.hash_;

  return *this;
}

auto UnorderedMapPool::operator[](const Key& key) -> Value& {
  size_type bucket_pos = bucket_hash(key);
  iterator it = find_hint(bucket_pos, key);

  if (it != end()) {
    return it->second;
  }
  
  if (size()+1 > threshold()) {
    rehash(bucket_count() * 2);
    bucket_pos = bucket_hash(key);
  }
  
  // not there -- insert a new one
  node* n = new node(std::make_pair(key, Value()));
  InternalInsert(buckets_, bucket_pos, n);

  ++size_;
  return n->value(); // n->data_.second
}

void UnorderedMapPool::rehash(size_type count) {
  node** tmp = new node*[count]();

  // iterate and rehash
  for (iterator it = begin(), last = end(); it != last; ++it) {
    // insert into tmp
    size_type hash = hash_function()(it->first) % count;
    push_front_bucket_node(tmp[hash], new node(*it));
  }
  buckets_ = tmp;
  bucket_count_ = count;
}

auto
UnorderedMapPool::insert(const value_type& value) -> std::pair<iterator, bool> {
  std::size_t bucket_pos = bucket_hash(value.first);
  iterator it = find_hint(bucket_pos, value.first);

  if (it != end()) {
    return std::make_pair(it, false);
  }

  it = InternalInsert(buckets_, bucket_pos, new node(value));
  ++size_;
  return std::make_pair(it, it != end());
}