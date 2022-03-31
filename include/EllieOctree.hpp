#ifndef OCTREE_H
#define OCTREE_H

#include <iostream>
#include <algorithm>
#include <iostream>
#include <utility>
#include <type_traits>
#include <concepts>

// Concepts
template<typename Tname>
concept copyable = std::is_copy_constructible_v<Tname>;

template<typename Tname>
concept copyableOnly = std::is_copy_constructible_v<Tname> && !std::is_move_constructible_v<Tname>;

template<typename Tname>
concept moveable = std::is_move_constructible_v<Tname>;

template<typename Tname>
concept moveableOnly = std::is_move_constructible_v<Tname> && !std::is_copy_constructible_v<Tname>;



struct vec3 {  // 3 vector struct for octree coordinates
	double x_, y_, z_;

	vec3 operator+(const vec3& V) const {  // +, addition
		vec3 v{ x_ + V.x_, y_ + V.y_, z_ + V.z_ };
		return v;
	}

	vec3 operator+(const double& S) const {  // +, addition
		vec3 v{ x_ + S, y_ + S, z_ + S };
		return v;
	}

	vec3 operator-(const vec3& V) const {  // -, subtraction
		vec3 v{ x_ - V.x_, y_ - V.y_, z_ - V.z_ };
		return v;
	}

	vec3 operator-(const double& S) const {  // -, subtraction
		vec3 v{ x_ - S, y_ - S, z_ - S };
		return v;
	}

	vec3 operator*(const double& S) const {  // *, multiplication, z1*z2
		vec3 v{ x_ * S, y_ * S, z_ * S };
		return v;
	}

	vec3 operator/(const double& S) const {  // /, division, z1/z2
		vec3 v{ x_ / S, y_ / S, z_ / S };
		return v;
	}

	friend vec3 operator+(const double S, const vec3& V) {
		vec3 v{ S + V.x_, S + V.y_, S + V.z_ };
		return v;
	}

	friend vec3 operator-(const double S, const vec3& V) {
		vec3 v{ S - V.x_, S - V.y_, S - V.z_ };
		return v;
	}

	friend vec3 operator*(const double S, const vec3& V) {
		vec3 v{ S * V.x_, S * V.y_, S * V.z_ };
		return v;
	}

	friend std::ostream& operator<<(std::ostream& os, const vec3& V) {  // overload << for particles
		os << "(" << V.x_ << ", " << V.y_ << ", " << V.z_ << ")";
		return os;
	}
};

template <typename Tname> class Node {  // Node structure
public:
	Tname* Objs_;  // an array of objects
	Node<Tname>* child_[8];  // pointers to child_ nodes
	Node<Tname>* parent_;  // parent_ pointer
	double x_, y_, z_;  // c.o.d. values
	double xMax_, xMin_, yMax_, yMin_, zMax_, zMin_;
	int num_;  // number of objects in node
	bool leaf_;  // bool of whether node is a leaf
	int depth_;  // depth of node in the tree

	// Rule of 5, with variety for objects that can't be copied/moved
	Node();
	template <copyableOnly T = Tname> Node(Tname* Objects, Node<Tname>* child0, Node<Tname>* child1, Node<Tname>* child2, Node<Tname>* child3,
		Node<Tname>* child4, Node<Tname>* child5, Node<Tname>* child6, Node<Tname>* child7, Node<Tname>* parent,
		double x, double y, double z, double xMax, double xMin, double yMax, double yMin, double zMax, double zMin,
		int num, int depth, bool leaf);
	template <moveable T = Tname> Node(Tname* Objects, Node<Tname>* child0, Node<Tname>* child1, Node<Tname>* child2, Node<Tname>* child3,
		Node<Tname>* child4, Node<Tname>* child5, Node<Tname>* child6, Node<Tname>* child7, Node<Tname>* parent,
		double x, double y, double z, double xMax, double xMin, double yMax, double yMin, double zMax, double zMin,
		int num, int depth, bool leaf);
	~Node();  // Destructor
	template <copyable T = Tname> Node(const Node<Tname>& node);  // Copy constructor
	template <moveableOnly T = Tname> Node(const Node<Tname>& node);  // Copy constructor
	template <copyable T = Tname> Node& operator = (const Node<Tname>& node);  // Copy Assignment operator
	template <moveableOnly T = Tname> Node& operator = (const Node<Tname>& node);  // Copy Assignment operator
	template <moveable T = Tname> Node(Node<Tname>&& node);  // Move constructor
	template <copyableOnly T = Tname> Node(Node<Tname>&& node);  // Move constructor
	template <moveable T = Tname> Node& operator = (Node<Tname>&& node); // Move assignment operator
	template <copyableOnly T = Tname> Node& operator = (Node<Tname>&& node); // Move assignment operator

	// setter/getter functions
	void setX(const double X) { x_ = X; }
	void setY(const double Y) { y_ = Y; }
	void setZ(const double Z) { z_ = Z; }
	double getXLength() { return xMax_ - xMin_; }
	double getYLength() { return yMax_ - yMin_; }
	double getZLength() { return zMax_ - zMin_; }
};

template <typename Tname>
class Octree {  // bound Octree container class
private:
	Node<Tname>* root_;  // pointer to the root Node
	int maxDepth_, maxLeaf_;  // maximum tree depth and leaf amount

public:
	typedef vec3(*objToCoord)(Tname& Obj);  // Function pointer
	objToCoord func;  // Function to convert an object to a coordinate

	// Rule of 5, with variety for objects that can't be copied/moved
	Octree();
	Octree(Tname* Objects, objToCoord funcToPlace, int objArrSize, int maxDepth, int maxLeafSize, double xMin = 0, double xMax = 0, double yMin = 0, double yMax = 0, double zMin = 0, double zMax = 0);  // Constructor
	~Octree();  // Destructor
	template <copyable T = Tname> Octree(const Octree<Tname>& O);  // Copy constructor
	template <moveableOnly T = Tname> Octree(const Octree<Tname>& O);  // Copy constructor
	template <copyable T = Tname> Octree& operator = (const Octree<Tname>& O);  // Copy Assignment operator
	template <moveableOnly T = Tname> Octree& operator = (const Octree<Tname>& O);  // Copy Assignment operator
	Octree(Octree<Tname>&& O);  // Move constructor
	Octree& operator = (Octree<Tname>&& O); // Move assignment operator

	// setter/getter functions
	void setFunc(objToCoord funcToPlace) { func = funcToPlace; }  // Sets the function to convert an object to a set of coords
	Node<Tname>* getRoot() const { return root_; };
	Node<Tname>* getNode(const int depth, const double X, const double Y, const double Z) const;  // Return Node at depth/position
	Tname* getNodeData(Node<Tname>* node) { return node->Objs_; }  // Return Node data
	double getMaxDepth() { return maxDepth_; }
	double getMaxLeaf() { return maxLeaf_; }
	int getDataSize(Node<Tname>* node, bool homeNode = true);  // Get data length at Node

	// Member functions
	Node<Tname>* findLeafNode(const double X, const double Y, const double Z) const;  // Returns a Leaf Node at position (x, y, z)
	template<copyableOnly T = Tname> void build(Node<Tname>* node, int depth = -1);  // Build tree
	template<moveable T = Tname> void build(Node<Tname>* node, int depth = -1);  // Build tree
	template<copyableOnly T = Tname> Tname* updateTree(Node<Tname>* node, Tname* lostObjs = nullptr, int oldRootSize = 0);
	template<moveable T = Tname> Tname* updateTree(Node<Tname>* node, Tname* lostObjs = nullptr, int oldRootSize = 0);
	template<copyableOnly T = Tname> void updateNode(Node<Tname>* node);
	template<moveable T = Tname> void updateNode(Node<Tname>* node);
	template<copyable T = Tname> Tname* copyTreeData(Node<Tname>* node, Tname* ObjArr = nullptr, bool homeNode = true) const;  // Return all data
	template<moveable T = Tname> Tname* moveTreeData(Node<Tname>* node, Tname* ObjArr = nullptr, bool homeNode = true);  // Return all data
	template<copyableOnly T = Tname> void addToTree(Tname Obj);  // Add an object
	template<moveable T = Tname> void addToTree(Tname Obj);  // Add an object
	template<copyableOnly T = Tname> void addToTree(Tname* ObjArr, int size);  // Add arrays of objects
	template<moveable T = Tname> void addToTree(Tname* ObjArr, int size);  // Add arrays of objects
};


// Node constructors/destructor
template <typename Tname> Node<Tname>::Node() {
	Objs_ = nullptr;
	for (int i = 0; i < 8; i++) {
		child_[i] = nullptr;
	}
	parent_ = nullptr;
	x_ = y_ = z_ = 0;
	num_ = 0;
	xMax_ = xMin_ = yMax_ = yMin_ = zMax_ = zMin_ = 0;
	depth_ = 0; leaf_ = false;
}

template <typename Tname>
template<copyableOnly>
Node<Tname>::Node(Tname* ObjArr, Node<Tname>* child0, Node<Tname>* child1, Node<Tname>* child2,
	Node<Tname>* child3, Node<Tname>* child4, Node<Tname>* child5,
	Node<Tname>* child6, Node<Tname>* child7, Node<Tname>* parent, double x,
	double y, double z, double xMax, double xMin, double yMax, double yMin,
	double zMax, double zMin, int num, int depth, bool leaf) {
	if (num > 0) {
		Objs_ = new Tname[num];
		for (int i = 0; i < num; i++) {
			Objs_[i] = ObjArr[i];
		}
	}
	else {
		Objs_ = nullptr;
	}
	child_[0] = child0; child_[1] = child1; child_[2] = child2; child_[3] = child3;
	child_[4] = child4; child_[5] = child5; child_[6] = child6; child_[7] = child7;
	parent_ = parent;
	x_ = x; y_ = y; z_ = z;
	xMax_ = xMax; xMin_ = xMin; yMax_ = yMax; yMin_ = yMin; zMax_ = zMax; zMin_ = zMin;
	num_ = num;
	depth_ = depth; leaf_ = leaf;
}

template <typename Tname>
template<moveable>
Node<Tname>::Node(Tname* ObjArr, Node<Tname>* child0, Node<Tname>* child1, Node<Tname>* child2,
	Node<Tname>* child3, Node<Tname>* child4, Node<Tname>* child5,
	Node<Tname>* child6, Node<Tname>* child7, Node<Tname>* parent, double x,
	double y, double z, double xMax, double xMin, double yMax, double yMin,
	double zMax, double zMin, int num, int depth, bool leaf) {
	if (num > 0) {
		Objs_ = new Tname[num];
		for (int i = 0; i < num; i++) {
			Objs_[i] = std::move(ObjArr[i]);
		}
	}
	else {
		Objs_ = nullptr;
	}
	child_[0] = child0; child_[1] = child1; child_[2] = child2; child_[3] = child3;
	child_[4] = child4; child_[5] = child5; child_[6] = child6; child_[7] = child7;
	parent_ = parent;
	x_ = x; y_ = y; z_ = z;
	xMax_ = xMax; xMin_ = xMin; yMax_ = yMax; yMin_ = yMin; zMax_ = zMax; zMin_ = zMin;
	num_ = num;
	depth_ = depth; leaf_ = leaf;
}

template <typename Tname> Node<Tname>::~Node() {
	if (Objs_) {  // If there is an array of objects
		delete[] Objs_;
	}
	for (int i = 0; i < 8; i++) {
		if (child_[i] != nullptr) { delete child_[i]; }  // iterate over children
	}
}

// Node copy constructors
template <typename Tname>
template <copyable>
Node<Tname>::Node(const Node<Tname>& node) {  // Node copy constructor
	Objs_ = new Tname[node.num_];  // Copy objects
	for (int i = 0; i < node.num_; i++) {
		Objs_[i] = node.Objs_[i];
	}
	for (int i = 0; i < 8; i++) {
		child_[i] = node.child_[i];
	}
	parent_ = node.parent_;
	x_ = node.x_; y_ = node.y_; z_ = node.z_;
	xMax_ = node.xMax_; xMin_ = node.xMin_; yMax_ = node.yMax_; yMin_ = node.yMin_; zMax_ = node.zMax_; zMin_ = node.zMin_;
	num_ = node.num_;
	depth_ = node.depth_;
	leaf_ = node.leaf_;
}

template <typename Tname>
template <moveableOnly>
Node<Tname>::Node(const Node<Tname>& node) {  // Node copy constructor
	std::cout << "Copying Node, moving Objects." << std::endl;
	Objs_ = new Tname[node.num_];  // Move objects
	for (int i = 0; i < node.num_; i++) {
		Objs_[i] = std::move(node.Objs_[i]);
	}
	delete[] node.Objs_;
	for (int i = 0; i < 8; i++) {
		child_[i] = node.child_[i];
	}
	parent_ = node.parent_;
	x_ = node.x_; y_ = node.y_; z_ = node.z_;
	xMax_ = node.xMax_; xMin_ = node.xMin_; yMax_ = node.yMax_; yMin_ = node.yMin_; zMax_ = node.zMax_; zMin_ = node.zMin_;
	num_ = node.num_;
	depth_ = node.depth_;
	leaf_ = node.leaf_;
}

// Node copy assign operators
template <typename Tname>
template <copyable>
Node<Tname>& Node<Tname>::operator = (const Node<Tname>& node) {  // Octree copy assignment operator
	if (&node == this) return *this;  // no self assignment
	delete[] Objs_;  // delete the object pointer array
	Objs_ = new Tname[node.num_];
	for (int i = 0; i < node.num_; i++) {
		Objs_[i] = node.Objs_[i];
	}
	delete[] child_;
	child_ = new Node[8];
	for (int i = 0; i < 8; i++) {
		child_[i] = node.child_[i];
	}
	delete[] parent_;
	parent_ = node.parent_;
	x_ = node.x_; y_ = node.y_; z_ = node.z_;
	xMax_ = node.xMax_; xMin_ = node.xMin_; yMax_ = node.yMax_; yMin_ = node.yMin_; zMax_ = node.zMax_; zMin_ = node.zMin_;
	num_ = node.num_;
	depth_ = node.depth_;
	leaf_ = node.leaf_;
	return *this;
}

template <typename Tname>
template <moveableOnly>
Node<Tname>& Node<Tname>::operator = (const Node<Tname>& node) {  // Node copy assignment operator
	std::cout << "Copying Node, moving Objects." << std::endl;
	if (&node == this) return *this;  // no self assignment
	delete[] Objs_;  // delete the object pointer array
	Objs_ = new Tname[node.num_];
	for (int i = 0; i < node.num_; i++) {
		Objs_[i] = std::move(node.Objs_[i]);
	}
	delete[] node.Objs_;
	delete[] child_;
	child_ = new Node[8];
	for (int i = 0; i < 8; i++) {
		child_[i] = node.child_[i];
	}
	delete[] parent_;
	parent_ = node.parent_;
	x_ = node.x_; y_ = node.y_; z_ = node.z_;
	xMax_ = node.xMax_; xMin_ = node.xMin_; yMax_ = node.yMax_; yMin_ = node.yMin_; zMax_ = node.zMax_; zMin_ = node.zMin_;
	num_ = node.num_;
	depth_ = node.depth_;
	leaf_ = node.leaf_;
	return *this;
}

// Node move constructors
template <typename Tname>
template <moveable>
Node<Tname>::Node(Node<Tname>&& node) {  // Node move constructor
	x_ = node.x_; y_ = node.y_; z_ = node.z_;
	xMax_ = node.xMax_; xMin_ = node.xMin_; yMax_ = node.yMax_; yMin_ = node.yMin_; zMax_ = node.zMax_; zMin_ = node.zMin_;
	num_ = node.num_;
	depth_ = node.depth_;
	leaf_ = node.leaf_;
	Objs_ = std::move(node.Objs_);
	parent_ = std::move(node.parent_);
	for (int i = 0; i < 8; i++) {
		child_[i] = std::move(node.child_[i]);
	}
	node.x_ = 0; node.y_ = 0; node.z_ = 0;
	node.xMax_ = 0; node.xMin_ = 0; node.yMax_ = 0; node.yMin_ = 0; node.zMax_ = 0; node.zMin_ = 0;
	node.num_ = 0;
	node.depth_ = 0;
	node.leaf_ = false;
	node.Objs_ = nullptr;
	node.parent_ = nullptr;
	for (int i = 0; i < 8; i++) {
		node.child_[i] = nullptr;
	}

}

template <typename Tname>
template <copyableOnly>
Node<Tname>::Node(Node<Tname>&& node) {  // Node move constructor
	std::cout << "Moving Node, copying Objects." << std::endl;
	x_ = node.x_; y_ = node.y_; z_ = node.z_;
	xMax_ = node.xMax_; xMin_ = node.xMin_; yMax_ = node.yMax_; yMin_ = node.yMin_; zMax_ = node.zMax_; zMin_ = node.zMin_;
	num_ = node.num_;
	depth_ = node.depth_;
	leaf_ = node.leaf_;
	Objs_ = node.Objs_;
	parent_ = std::move(node.parent_);
	for (int i = 0; i < 8; i++) {
		child_[i] = std::move(node.child_[i]);
	}
	node.x_ = 0; node.y_ = 0; node.z_ = 0;
	node.xMax_ = 0; node.xMin_ = 0; node.yMax_ = 0; node.yMin_ = 0; node.zMax_ = 0; node.zMin_ = 0;
	node.num_ = 0;
	node.depth_ = 0;
	node.leaf_ = false;
	node.Objs_ = nullptr;
	node.parent_ = nullptr;
	for (int i = 0; i < 8; i++) {
		node.child_[i] = nullptr;
	}

}

// Node move assign operators
template <typename Tname>
template <moveable>
Node<Tname>& Node<Tname>::operator = (Node<Tname>&& node) {  // Node move assignment operator
	std::swap(x_, node.x_); std::swap(y_, node.y_); std::swap(z_, node.z_);
	std::swap(xMax_, node.xMax_); std::swap(yMax_, node.yMax_); std::swap(zMax_, node.zMax_);
	std::swap(xMin_, node.xMin_); std::swap(yMin_, node.yMin_); std::swap(zMin_, node.zMin_);
	std::swap(num_, node.num_);
	std::swap(depth_, node.depth_);
	std::swap(leaf_, node.leaf_);
	std::swap(Objs_, node.Objs_);
	std::swap(parent_, node.parent_);
	for (int i = 0; i < 8; i++) {
		std::swap(child_[i], node.child_[i]);
	}
	return *this;
}

template <typename Tname>
template <copyableOnly>
Node<Tname>& Node<Tname>::operator = (Node<Tname>&& node) {  // Node move assignment operator
	std::cout << "Moving Node, copying Objects." << std::endl;
	std::swap(x_, node.x_); std::swap(y_, node.y_); std::swap(z_, node.z_);
	std::swap(xMax_, node.xMax_); std::swap(yMax_, node.yMax_); std::swap(zMax_, node.zMax_);
	std::swap(xMin_, node.xMin_); std::swap(yMin_, node.yMin_); std::swap(zMin_, node.zMin_);
	std::swap(num_, node.num_);
	std::swap(depth_, node.depth_);
	std::swap(leaf_, node.leaf_);
	Objs_ = node.Objs_; delete[] node.Objs_; node.Objs_ = nullptr;
	std::swap(parent_, node.parent_);
	for (int i = 0; i < 8; i++) {
		std::swap(child_[i], node.child_[i]);
	}
	return *this;
}



// Octree constructors/destructor
template <typename Tname> Octree<Tname>::Octree() {
	root_ = new Node<Tname>();
	maxDepth_ = 0;
	maxLeaf_ = 0;
}

template <typename Tname> Octree<Tname>::Octree(Tname* Objects, objToCoord funcToPlace, int objArrSize, int maxDepth, int maxLeafSize, double xMin, double xMax, double yMin, double yMax, double zMin, double zMax) {
	setFunc(funcToPlace);
	maxDepth_ = maxDepth;
	maxLeaf_ = maxLeafSize;
	root_ = new Node<Tname>(Objects, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0, 0, 0, xMax, xMin, yMax, yMin, zMax, zMin, objArrSize, -1, false);  // make root
//	*root_ = { nullptr,{ nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr }, nullptr, 0, 0, 0, 0, 0, 0, 0, 0, 0, objArrSize, 0 };  // initialise Node
//	root_->Objs_ = Objects;  // assign root object array
	build(root_);  // build root children
}

template <typename Tname> Octree<Tname>::~Octree() {  // Default destructor
	if (root_) {
		delete root_;  // calls Node destructor
	}
}

// Octree copy constructors
template <typename Tname>
template <copyable>
Octree<Tname>::Octree(const Octree<Tname>& O) {  // Octree copy constructor
	root_ = 0; root_ = O.getRoot();  // Copy root
	maxDepth_ = O.maxDepth_; maxLeaf_ = O.maxLeaf_;
	func = O.func;  // set function
	Tname* ObjArr = O.copyTreeData(root_);  // get object array (through copy)
	root_->Objs_ = ObjArr;
	build(root_);  // build tree
}

template <typename Tname>
template <moveableOnly>
Octree<Tname>::Octree(const Octree<Tname>& O) {  // Octree copy constructor
	std::cout << "Copying Octree, moving Objects." << std::endl;
	root_ = 0; root_ = O.getRoot();  // Copy root
	maxDepth_ = O.maxDepth_; maxLeaf_ = O.maxLeaf_;
	func = O.func;  // set function
	Tname* ObjArr = O.moveTreeData(root_);  // get object array (through move)
	root_->Objs_ = ObjArr;
	build(root_);  // build tree
}

// Octree copy assign operators
template <typename Tname>
template <copyable>
Octree<Tname>& Octree<Tname>::operator = (const Octree<Tname>& O) {  // Octree copy assignment operator
	if (&O == this) return *this;  // no self assignment
	delete root_;
	root_ = nullptr;
	maxDepth_ = 0; maxLeaf_ = 0; // delete this object�s array
	Tname* ObjArr = O.copyTreeData(O.getRoot());
	root_ = new Node<Tname>(ObjArr, nullptr, nullptr, nullptr,
		nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		O.getRoot()->num_, 0, 0, false);
	delete[] ObjArr;
	maxDepth_ = O.maxDepth_; maxLeaf_ = O.maxLeaf_;
	func = O.func;
	build(root_);
	std::cout << "Copying Octree." << std::endl;
	return *this;
}

template <typename Tname>
template <moveableOnly>
Octree<Tname>& Octree<Tname>::operator = (const Octree<Tname>& O) {  // Octree copy assignment operator
	if (&O == this) return *this;  // no self assignment
	delete root_;
	root_ = nullptr;
	maxDepth_ = 0; maxLeaf_ = 0; // delete this object�s array
	Tname* ObjArr = O.moveTreeData(O.getRoot());
	root_ = new Node<Tname>(ObjArr, nullptr, nullptr, nullptr,
		nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		O.getRoot()->num_, 0, 0, false);
	delete[] ObjArr;
	maxDepth_ = O.maxDepth_; maxLeaf_ = O.maxLeaf_;
	func = O.func;
	build(root_);
	std::cout << "Copying Octree, moving Objects." << std::endl;
	return *this;
}

// Octree move constructor & assign operator
template <typename Tname> Octree<Tname>::Octree(Octree<Tname>&& O) {  // Octree move constructor
	std::cout << "Move constructor for Octree called" << std::endl;
	maxDepth_ = O.maxDepth_; maxLeaf_ = O.maxLeaf_;
	root_ = O.root_;
	O.maxDepth_ = 0; O.maxLeaf_ = 0;
	O.root_ = 0;
}

template <typename Tname> Octree<Tname>& Octree<Tname>::operator = (Octree<Tname>&& O) {  // Octree move assignment operator
	std::cout << "Octree move assignment" << std::endl;
	std::swap(maxDepth_, O.maxDepth_); std::swap(maxLeaf_, O.maxLeaf_);
	std::swap(root_, O.root_);
	return *this;
}

// Octree getter functions
template <typename Tname> Node<Tname>* Octree<Tname>::getNode(const int depth, const double X, const double Y, const double Z) const {
	int curdepth = 0;
	Node<Tname>* node = root_;
	if (depth < maxDepth_) {
		while ((curdepth != depth) && (node->num_ > maxLeaf_)) {  // while not at required depth
			double xCent = (node->xMax_ + node->xMin_) / 2;
			double yCent = (node->yMax_ + node->yMin_) / 2;
			double zCent = (node->zMax_ + node->zMin_) / 2;
			bool front = X <= xCent;  // Determine node to travel to
			bool left = Y <= yCent;
			bool bottom = Z <= zCent;
			if (left && front && bottom) {  // child_[0]
				node = node->child_[0];
			}
			else if (!left && front && bottom) {  // 1
				node = node->child_[1];
			}
			else if (left && !front && bottom) {  // 2
				node = node->child_[2];
			}
			else if (!left && !front && bottom) {  // 3
				node = node->child_[3];
			}
			else if (left && front && !bottom) {  // 4
				node = node->child_[4];
			}
			else if (!left && front && !bottom) {  // 5
				node = node->child_[5];
			}
			else if (left && !front && !bottom) {  // 6
				node = node->child_[6];
			}
			else {  // 7
				node = node->child_[7];
			}
			curdepth++;
		}
		if (curdepth != depth) {
			std::cout << "You are at a leaf Node. You cannot go any deeper." << std::endl;
		}
	}
	else {
		std::cout << "Depth is too low. Please keep depth to " << maxDepth_ << " or less." << std::endl;
	}
	return node;
}

template <typename Tname> int Octree<Tname>::getDataSize(Node<Tname>* node, bool homeNode) {
	static int size;  // return the number of objects in a node's children
	if (homeNode) {
		homeNode = false;
		size = 0;
	}
	for (int i = 0; i < 8; i++) {  // Iterating over the child_ nodes
		if (node->child_[i]->num_ > maxLeaf_ && node->child_[i]->depth_ != maxDepth_) {
			getDataSize(node->child_[i], homeNode);  // repeat process with child_ nodes
		}
		else {
			size = size + node->child_[i]->num_;
		}
	}
	return size;
}


template <typename Tname> Node<Tname>* Octree<Tname>::findLeafNode(const double X, const double Y, const double Z) const {
	Node<Tname>* node = root_;
	if (X > node->xMax_ || X < node->xMin_ || Y > node->yMax_ || Y < node->yMin_ || Z > node->zMax_ || Z < node->zMin_) {
		return nullptr; // Coordinates are out of bounds of the octree, return nullptr
	}
	while (node->leaf_ == false) {  // while not at a leaf node
		if (node->num_ == 0) { break; }  // At an empty node, return the empty node
		double xCent = (node->xMax_ + node->xMin_) / 2;
		double yCent = (node->yMax_ + node->yMin_) / 2;
		double zCent = (node->zMax_ + node->zMin_) / 2;
		bool front = X <= xCent;  // Determine node to travel to
		bool left = Y <= yCent;
		bool bottom = Z <= zCent;
		if (left) {
			if (front) {
				if (bottom) {
					node = node->child_[0];
				}
				else {
					node = node->child_[4];
				}
			}
			else {
				if (bottom) {
					node = node->child_[2];
				}
				else {
					node = node->child_[6];
				}
			}
		}
		else {
			if (front) {
				if (bottom) {
					node = node->child_[1];
				}
				else {
					node = node->child_[5];
				}
			}
			else {
				if (bottom) {
					node = node->child_[3];
				}
				else {
					node = node->child_[7];
				}
			}
		}
	}
	return node;
}

// Octee build functions
template <typename Tname>
template<copyableOnly>
void Octree<Tname>::build(Node<Tname>* node, int depth) {
	depth++;
	vec3* V = new vec3[node->num_];
	for (int i = 0; i < node->num_; i++) {  // Calculate c.o.d. of Node by producing array of vec3=length of object array
		V[i] = func(node->Objs_[i]);  // using function provided get obj coords
		if (node == root_) {
			if (i == 0) {  // First iteration, set limits
				if (root_->xMin_ == root_->xMax_ && root_->yMin_ == root_->yMax_ && root_->zMin_ == root_->zMax_) {  // If root limits aren't preset
					root_->xMax_ = V[i].x_; root_->xMin_ = V[i].x_;
					root_->yMax_ = V[i].y_; root_->yMin_ = V[i].y_;
					root_->zMax_ = V[i].z_; root_->zMin_ = V[i].z_;
				}
			}  // Test if objects exceed current root bounds
			root_->xMax_ = (V[i].x_ > root_->xMax_) ? V[i].x_ : root_->xMax_;
			root_->xMin_ = (V[i].x_ < root_->xMin_) ? V[i].x_ : root_->xMin_;
			root_->yMax_ = (V[i].y_ > root_->yMax_) ? V[i].y_ : root_->yMax_;
			root_->yMin_ = (V[i].y_ < root_->yMin_) ? V[i].y_ : root_->yMin_;
			root_->zMax_ = (V[i].z_ > root_->zMax_) ? V[i].z_ : root_->zMax_;
			root_->zMin_ = (V[i].z_ < root_->zMin_) ? V[i].z_ : root_->zMin_;
		}
		node->x_ += V[i].x_;
		node->y_ += V[i].y_;
		node->z_ += V[i].z_;
	}
	if (node->num_ > 0) {
		node->setX(node->x_ / node->num_); node->setY(node->y_ / node->num_); node->setZ(node->z_ / node->num_);
	}
	double xCent = (node->xMax_ + node->xMin_) / 2;
	double yCent = (node->yMax_ + node->yMin_) / 2;
	double zCent = (node->zMax_ + node->zMin_) / 2;

	int Total = node->num_;
	Tname* Ob0, * Ob1, * Ob2, * Ob3, * Ob4, * Ob5, * Ob6, * Ob7;
	int BLF = 0, BRF = 0, BLB = 0, BRB = 0, TLF = 0, TRF = 0, TLB = 0, TRB = 0;
	for (int i = 0; i < Total; i++) {  // Get sizes of each box
		bool front = V[i].x_ <= xCent;
		bool left = V[i].y_ <= yCent;
		bool bottom = V[i].z_ <= zCent;
		if (left && front && bottom) {
			BLF++;
		}
		else if (!left && front && bottom) {
			BRF++;
		}
		else if (left && !front && bottom) {
			BLB++;
		}
		else if (!left && !front && bottom) {
			BRB++;
		}
		else if (left && front && !bottom) {
			TLF++;
		}
		else if (!left && front && !bottom) {
			TRF++;
		}
		else if (left && !front && !bottom) {
			TLB++;
		}
		else {
			TRB++;
		}
	}  // replicate numbers and initialise arrays
	const int num0 = BLF, num1 = BRF, num2 = BLB, num3 = BRB, num4 = TLF, num5 = TRF, num6 = TLB, num7 = TRB;
	Ob0 = new Tname[BLF]; Ob1 = new Tname[BRF]; Ob2 = new Tname[BLB]; Ob3 = new Tname[BRB];
	Ob4 = new Tname[TLF]; Ob5 = new Tname[TRF]; Ob6 = new Tname[TLB]; Ob7 = new Tname[TRB];
	for (int i = 0; i < Total; i++) {  // Copy over relevant particles
		bool front = V[i].x_ <= xCent;
		bool left = V[i].y_ <= yCent;
		bool bottom = V[i].z_ <= zCent;
		if (left && front && bottom) {  // particles sent to correct child_ node
			Ob0[BLF - 1] = node->Objs_[i];
			BLF--;
		}
		else if (!left && front && bottom) {
			Ob1[BRF - 1] = node->Objs_[i];
			BRF--;
		}
		else if (left && !front && bottom) {
			Ob2[BLB - 1] = node->Objs_[i];
			BLB--;
		}
		else if (!left && !front && bottom) {
			Ob3[BRB - 1] = node->Objs_[i];
			BRB--;
		}
		else if (left && front && !bottom) {
			Ob4[TLF - 1] = node->Objs_[i];
			TLF--;
		}
		else if (!left && front && !bottom) {
			Ob5[TRF - 1] = node->Objs_[i];
			TRF--;
		}
		else if (left && !front && !bottom) {
			Ob6[TLB - 1] = node->Objs_[i];
			TLB--;
		}
		else {
			Ob7[TRB - 1] = node->Objs_[i];
			TRB--;
		}
	}
	// Create and initialise child_ nodes with pointers to the parent_ node and relevant object arrays
	// Node(Objects, child0, child1, child2, child3, child4, child5, child6, child7, parent, x, y, z, xMax, xMin, yMax, yMin, zMax, zMin, num, depth, leaf);
	node->child_[0] = new Node<Tname>(Ob0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0, 0, 0, xCent, node->xMin_, yCent, node->yMin_, zCent, node->zMin_, num0, depth, false);
	node->child_[1] = new Node<Tname>(Ob1, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0, 0, 0, xCent, node->xMin_, node->yMax_, yCent, zCent, node->zMin_, num1, depth, false);
	node->child_[2] = new Node<Tname>(Ob2, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0, 0, 0, node->xMax_, xCent, yCent, node->yMin_, zCent, node->zMin_, num2, depth, false);
	node->child_[3] = new Node<Tname>(Ob3, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0, 0, 0, node->xMax_, xCent, node->yMax_, yCent, zCent, node->zMin_, num3, depth, false);
	node->child_[4] = new Node<Tname>(Ob4, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0, 0, 0, xCent, node->xMin_, yCent, node->yMin_, node->zMax_, zCent, num4, depth, false);
	node->child_[5] = new Node<Tname>(Ob5, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0, 0, 0, xCent, node->xMin_, node->yMax_, yCent, node->zMax_, zCent, num5, depth, false);
	node->child_[6] = new Node<Tname>(Ob6, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0, 0, 0, node->xMax_, xCent, yCent, node->yMin_, node->zMax_, zCent, num6, depth, false);
	node->child_[7] = new Node<Tname>(Ob7, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0, 0, 0, node->xMax_, xCent, node->yMax_, yCent, node->zMax_, zCent, num7, depth, false);
	node->child_[0]->parent_ = node;
	node->child_[1]->parent_ = node;
	node->child_[2]->parent_ = node;
	node->child_[3]->parent_ = node;
	node->child_[4]->parent_ = node;
	node->child_[5]->parent_ = node;
	node->child_[6]->parent_ = node;
	node->child_[7]->parent_ = node;
	delete[] V;  // delete array of obj coords
	delete[] Ob0; delete[] Ob1; delete[] Ob2; delete[] Ob3; delete[] Ob4; delete[] Ob5; delete[] Ob6; delete[] Ob7;
	for (int i = 0; i < 8; i++) {
		if ((node->child_[i]->num_ > maxLeaf_) && (depth != maxDepth_)) {  // If there is a significant number of particles in box and not at max depth
			build(node->child_[i], depth);  // repeat process with child_ nodes
		}
		else {  // child_ node is a leaf node or an empty node
			if (node->child_[i]->num_ == 0) {  // empty case
				node->child_[i]->Objs_ = nullptr;
				continue;
			}  // leaf case
			vec3* V = new vec3[node->child_[i]->num_];  // get child_ c.o.d. coords
			std::pair<vec3, double> result;
			for (int j = 0; j < node->child_[i]->num_; j++) {  // Calculate c.o.d. of node by producing array of vec3=length of object array
				V[j] = func(node->child_[i]->Objs_[j]);
				node->child_[i]->x_ += V[j].x_;
				node->child_[i]->y_ += V[j].y_;
				node->child_[i]->z_ += V[j].z_;
			}
			node->child_[i]->leaf_ = true;
			delete[] V;
			if (node->child_[i]->num_ > 0) {  // set c.o.d.
				node->child_[i]->setX(node->child_[i]->x_ / node->child_[i]->num_);
				node->child_[i]->setY(node->child_[i]->y_ / node->child_[i]->num_);
				node->child_[i]->setZ(node->child_[i]->z_ / node->child_[i]->num_);
			}
		}
	}

	if (node->num_ > maxLeaf_) {  // If there is a significant number of particles in box
		if (depth != maxDepth_) {  // and depth is not child_ node depth (and hence definitely not a leaf node)
			if (node->Objs_) { delete[] node->Objs_; }  // delete node objects, freeing up memory
			node->Objs_ = nullptr;
		}
	}
}

template <typename Tname>
template<moveable>
void Octree<Tname>::build(Node<Tname>* node, int depth) {
	depth++;
	vec3* V = new vec3[node->num_];
	for (int i = 0; i < node->num_; i++) {  // Calculate c.o.d. of Node by producing array of vec3=length of object array
		V[i] = func(node->Objs_[i]);  // using function provided get obj coords
		if (node == root_) {
			if (i == 0) {  // First iteration, set limits
				if (root_->xMin_ == root_->xMax_ && root_->yMin_ == root_->yMax_ && root_->zMin_ == root_->zMax_) {  // If root limits aren't preset
					root_->xMax_ = V[i].x_; root_->xMin_ = V[i].x_;
					root_->yMax_ = V[i].y_; root_->yMin_ = V[i].y_;
					root_->zMax_ = V[i].z_; root_->zMin_ = V[i].z_;
				}
			}  // Test if objects exceed current root bounds
			root_->xMax_ = (V[i].x_ > root_->xMax_) ? V[i].x_ : root_->xMax_;
			root_->xMin_ = (V[i].x_ < root_->xMin_) ? V[i].x_ : root_->xMin_;
			root_->yMax_ = (V[i].y_ > root_->yMax_) ? V[i].y_ : root_->yMax_;
			root_->yMin_ = (V[i].y_ < root_->yMin_) ? V[i].y_ : root_->yMin_;
			root_->zMax_ = (V[i].z_ > root_->zMax_) ? V[i].z_ : root_->zMax_;
			root_->zMin_ = (V[i].z_ < root_->zMin_) ? V[i].z_ : root_->zMin_;
		}
		node->x_ += V[i].x_;
		node->y_ += V[i].y_;
		node->z_ += V[i].z_;
	}
	if (node->num_ > 0) {
		node->setX(node->x_ / node->num_); node->setY(node->y_ / node->num_); node->setZ(node->z_ / node->num_);
	}
	double xCent = (node->xMax_ + node->xMin_) / 2;
	double yCent = (node->yMax_ + node->yMin_) / 2;
	double zCent = (node->zMax_ + node->zMin_) / 2;

	int Total = node->num_;
	Tname* Ob0, * Ob1, * Ob2, * Ob3, * Ob4, * Ob5, * Ob6, * Ob7;
	int BLF = 0, BRF = 0, BLB = 0, BRB = 0, TLF = 0, TRF = 0, TLB = 0, TRB = 0;
	for (int i = 0; i < Total; i++) {  // Get sizes of each box
		bool front = V[i].x_ <= xCent;
		bool left = V[i].y_ <= yCent;
		bool bottom = V[i].z_ <= zCent;
		if (left && front && bottom) {
			BLF++;
		}
		else if (!left && front && bottom) {
			BRF++;
		}
		else if (left && !front && bottom) {
			BLB++;
		}
		else if (!left && !front && bottom) {
			BRB++;
		}
		else if (left && front && !bottom) {
			TLF++;
		}
		else if (!left && front && !bottom) {
			TRF++;
		}
		else if (left && !front && !bottom) {
			TLB++;
		}
		else {
			TRB++;
		}
	}  // replicate numbers and initialise arrays
	const int num0 = BLF, num1 = BRF, num2 = BLB, num3 = BRB, num4 = TLF, num5 = TRF, num6 = TLB, num7 = TRB;
	Ob0 = new Tname[BLF]; Ob1 = new Tname[BRF]; Ob2 = new Tname[BLB]; Ob3 = new Tname[BRB];
	Ob4 = new Tname[TLF]; Ob5 = new Tname[TRF]; Ob6 = new Tname[TLB]; Ob7 = new Tname[TRB];
	for (int i = 0; i < Total; i++) {  // Copy over relevant particles
		bool front = V[i].x_ <= xCent;
		bool left = V[i].y_ <= yCent;
		bool bottom = V[i].z_ <= zCent;
		if (left && front && bottom) {  // particles sent to correct child_ node
			Ob0[BLF - 1] = std::move(node->Objs_[i]);
			BLF--;
		}
		else if (!left && front && bottom) {
			Ob1[BRF - 1] = std::move(node->Objs_[i]);
			BRF--;
		}
		else if (left && !front && bottom) {
			Ob2[BLB - 1] = std::move(node->Objs_[i]);
			BLB--;
		}
		else if (!left && !front && bottom) {
			Ob3[BRB - 1] = std::move(node->Objs_[i]);
			BRB--;
		}
		else if (left && front && !bottom) {
			Ob4[TLF - 1] = std::move(node->Objs_[i]);
			TLF--;
		}
		else if (!left && front && !bottom) {
			Ob5[TRF - 1] = std::move(node->Objs_[i]);
			TRF--;
		}
		else if (left && !front && !bottom) {
			Ob6[TLB - 1] = std::move(node->Objs_[i]);
			TLB--;
		}
		else {
			Ob7[TRB - 1] = std::move(node->Objs_[i]);
			TRB--;
		}
	}
	// Create and initialise child_ nodes with pointers to the parent_ node and relevant object arrays
	// Node(Objects, child0, child1, child2, child3, child4, child5, child6, child7, parent, x, y, z, xMax, xMin, yMax, yMin, zMax, zMin, num, depth, leaf);
	node->child_[0] = new Node<Tname>(Ob0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0, 0, 0, xCent, node->xMin_, yCent, node->yMin_, zCent, node->zMin_, num0, depth, false);
	node->child_[1] = new Node<Tname>(Ob1, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0, 0, 0, xCent, node->xMin_, node->yMax_, yCent, zCent, node->zMin_, num1, depth, false);
	node->child_[2] = new Node<Tname>(Ob2, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0, 0, 0, node->xMax_, xCent, yCent, node->yMin_, zCent, node->zMin_, num2, depth, false);
	node->child_[3] = new Node<Tname>(Ob3, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0, 0, 0, node->xMax_, xCent, node->yMax_, yCent, zCent, node->zMin_, num3, depth, false);
	node->child_[4] = new Node<Tname>(Ob4, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0, 0, 0, xCent, node->xMin_, yCent, node->yMin_, node->zMax_, zCent, num4, depth, false);
	node->child_[5] = new Node<Tname>(Ob5, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0, 0, 0, xCent, node->xMin_, node->yMax_, yCent, node->zMax_, zCent, num5, depth, false);
	node->child_[6] = new Node<Tname>(Ob6, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0, 0, 0, node->xMax_, xCent, yCent, node->yMin_, node->zMax_, zCent, num6, depth, false);
	node->child_[7] = new Node<Tname>(Ob7, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0, 0, 0, node->xMax_, xCent, node->yMax_, yCent, node->zMax_, zCent, num7, depth, false);
	node->child_[0]->parent_ = node;
	node->child_[1]->parent_ = node;
	node->child_[2]->parent_ = node;
	node->child_[3]->parent_ = node;
	node->child_[4]->parent_ = node;
	node->child_[5]->parent_ = node;
	node->child_[6]->parent_ = node;
	node->child_[7]->parent_ = node;
	delete[] V;  // delete array of obj coords
	delete[] Ob0; delete[] Ob1; delete[] Ob2; delete[] Ob3; delete[] Ob4; delete[] Ob5; delete[] Ob6; delete[] Ob7;
	for (int i = 0; i < 8; i++) {
		if ((node->child_[i]->num_ > maxLeaf_) && (depth != maxDepth_)) {  // If there is a significant number of particles in box and not at max depth
			build(node->child_[i], depth);  // repeat process with child_ nodes
		}
		else {  // child_ node is a leaf node or an empty node
			if (node->child_[i]->num_ == 0) {  // empty case
				node->child_[i]->Objs_ = nullptr;
				continue;
			}  // leaf case
			vec3* V = new vec3[node->child_[i]->num_];  // get child_ c.o.d. coords
			std::pair<vec3, double> result;
			for (int j = 0; j < node->child_[i]->num_; j++) {  // Calculate c.o.d. of node by producing array of vec3=length of object array
				V[j] = func(node->child_[i]->Objs_[j]);
				node->child_[i]->x_ += V[j].x_;
				node->child_[i]->y_ += V[j].y_;
				node->child_[i]->z_ += V[j].z_;
			}
			node->child_[i]->leaf_ = true;
			delete[] V;
			if (node->child_[i]->num_ > 0) {  // set c.o.d.
				node->child_[i]->setX(node->child_[i]->x_ / node->child_[i]->num_);
				node->child_[i]->setY(node->child_[i]->y_ / node->child_[i]->num_);
				node->child_[i]->setZ(node->child_[i]->z_ / node->child_[i]->num_);
			}
		}
	}

	if (node->num_ > maxLeaf_) {  // If there is a significant number of particles in box
		if (depth != maxDepth_) {  // and depth is not child_ node depth (and hence definitely not a leaf node)
			if (node->Objs_) { delete[] node->Objs_; }  // delete node objects, freeing up memory
			node->Objs_ = nullptr;
		}
	}
}

// Octree update functions
template <typename Tname>
template <copyableOnly>
Tname* Octree<Tname>::updateTree(Node<Tname>* node, Tname* lostObjs, int oldRootSize) {
	// Iterate over nodes
	if (oldRootSize == 0) { oldRootSize = root_->num_; }
	int lostObjsSize = oldRootSize - root_->num_;
	for (int i = 0; i < 8; i++) {
		if (node->child_[i] && node->child_[i]->num_ > 0) {
			if (node->child_[i]->leaf_) {	// child[i] is leaf node
				// For each object in a leaf node, check coordinates and find the new leaf node
				int oldNum = node->child_[i]->num_; int newNum = oldNum;
				bool* leftObjIndx = new bool[oldNum]();
				for (int j = 0; j < node->child_[i]->num_; j++) {
					leftObjIndx[j] = false;
					vec3 coords = func(node->child_[i]->Objs_[j]);  // using function provided to get obj coords
					Node<Tname>* destNode = findLeafNode(coords.x_, coords.y_, coords.z_);
					if (!destNode) {  // If object leaves octree, add to lostObjs
						if (lostObjsSize != 0) {
							lostObjsSize++;
							Tname* newLostObjs = new Tname[lostObjsSize]();
							for (int k = 0; k < lostObjsSize - 1; k++) {
								newLostObjs[k] = lostObjs[k];
							}
							newLostObjs[lostObjsSize - 1] = node->child_[i]->Objs_[j];
							delete[] lostObjs;
							lostObjs = newLostObjs;
						}
						else {
							lostObjsSize++;
							lostObjs = new Tname[1];
							lostObjs[0] = node->child_[i]->Objs_[j];
						}
						// Remove object from leaf node
						leftObjIndx[j] = true;
						newNum--;
					}
					else if (destNode != node->child_[i]) {  // If destination node!=origin node, move object to new node
						destNode->num_++;
						Tname* destObjs = new Tname[destNode->num_]();
						for (int k = 0; k < destNode->num_ - 1; k++) {
							destObjs[k] = destNode->Objs_[k];
						}
						// Move object to leaf node, adjust array sizes appropriately
						destObjs[destNode->num_ - 1] = node->child_[i]->Objs_[j];
						if (destNode->Objs_) {
							delete[] destNode->Objs_;
						}
						destNode->Objs_ = destObjs;
						updateNode(destNode);  // Update destination node, building more leaves if needed, update extended family nodes
						// Remove object from leaf node
						leftObjIndx[j] = true;
						newNum--;
					}
				}
				if (oldNum != newNum) {  // If object(s) have left child node
					if (newNum == 0) {  // child node is empty
						node->child_[i]->num_ = newNum;
						updateNode(node->child_[i]);  // Update node and its parents
					}
					else {  // child node is still a leaf
						Tname* originObjs = new Tname[newNum];
						int counter = 0;
						for (int j = 0; j < oldNum; j++) {  // Remove objects that have left child node
							if (leftObjIndx[j] == false) {
								originObjs[counter] = node->child_[i]->Objs_[j];
								counter++;
							}
						}
						node->child_[i]->num_ = newNum;
						delete[] node->child_[i]->Objs_;
						node->child_[i]->Objs_ = originObjs;
						updateNode(node->child_[i]);  // Update child's parent nodes
					}
				}
				lostObjsSize = oldRootSize - root_->num_;
				delete[] leftObjIndx;
			}
			else {  // repeat for child's child nodes
				lostObjs = updateTree(node->child_[i], lostObjs, oldRootSize);
			}
		}
	}
	return lostObjs;
}

template <typename Tname>
template <moveable>
Tname* Octree<Tname>::updateTree(Node<Tname>* node, Tname* lostObjs, int oldRootSize) {
	// Iterate over nodes
	if (oldRootSize == 0) { oldRootSize = root_->num_; }
	int lostObjsSize = oldRootSize - root_->num_;
	for (int i = 0; i < 8; i++) {
		if (node->child_[i] && node->child_[i]->num_ > 0) {
			if (node->child_[i]->leaf_) {	// child[i] is leaf node
				// For each object in a leaf node, check coordinates and find the new leaf node
				int oldNum = node->child_[i]->num_; int newNum = oldNum;
				bool* leftObjIndx = new bool[oldNum]();
				for (int j = 0; j < node->child_[i]->num_; j++) {
					leftObjIndx[j] = false;
					vec3 coords = func(node->child_[i]->Objs_[j]);  // using function provided to get obj coords
					Node<Tname>* destNode = findLeafNode(coords.x_, coords.y_, coords.z_);
					if (!destNode) {  // If object leaves octree, add to lostObjs
						if (lostObjsSize != 0) {
							lostObjsSize++;
							Tname* newLostObjs = new Tname[lostObjsSize];
							for (int k = 0; k < lostObjsSize - 1; k++) {
								newLostObjs[k] = std::move(lostObjs[k]);
							}
							newLostObjs[lostObjsSize - 1] = std::move(node->child_[i]->Objs_[j]);
							delete[] lostObjs;
							lostObjs = newLostObjs;
						}
						else {
							lostObjsSize++;
							lostObjs = new Tname[1];
							lostObjs[0] = std::move(node->child_[i]->Objs_[j]);
						}
						// Remove object from leaf node
						leftObjIndx[j] = true;
						newNum--;
					}
					else if (destNode != node->child_[i]) {  // If destination node!=origin node, move object to new node
						destNode->num_++;
						Tname* destObjs = new Tname[destNode->num_];
						for (int k = 0; k < destNode->num_ - 1; k++) {
							destObjs[k] = std::move(destNode->Objs_[k]);
						}
						// Move object to leaf node, adjust array sizes appropriately
						destObjs[destNode->num_ - 1] = std::move(node->child_[i]->Objs_[j]);
						if (destNode->Objs_) {
							delete[] destNode->Objs_;
						}
						destNode->Objs_ = destObjs;
						updateNode(destNode);  // Update destination node, building more leaves if needed, update extended family nodes
						// Remove object from leaf node
						leftObjIndx[j] = true;
						newNum--;
					}
				}
				if (oldNum != newNum) {  // If object(s) have left child node
					if (newNum == 0) {  // child node is empty
						node->child_[i]->num_ = newNum;
						updateNode(node->child_[i]);  // Update node and its parents
					}
					else {  // child node is still a leaf
						Tname* originObjs = new Tname[newNum];
						int counter = 0;
						for (int j = 0; j < oldNum; j++) {  // Remove objects that have left child node
							if (leftObjIndx[j] == false) {
								originObjs[counter] = std::move(node->child_[i]->Objs_[j]);
								counter++;
							}
						}
						node->child_[i]->num_ = newNum;
						delete[] node->child_[i]->Objs_;
						node->child_[i]->Objs_ = originObjs;
						updateNode(node->child_[i]);  // Update child's parent nodes
					}
				}
				lostObjsSize = oldRootSize - root_->num_;
				delete[] leftObjIndx;
			}
			else {  // repeat for child's child nodes
				lostObjs = updateTree(node->child_[i], lostObjs, oldRootSize);
			}
		}
	}
	return lostObjs;
}

template <typename Tname>
template <copyableOnly>
void Octree<Tname>::updateNode(Node<Tname>* node) {
	// Updates node statistics with information from object list
	if (node->leaf_) {  // If node was previously a leaf node
						// Node has gained enough particles and is no longer a leaf
						// Node has different number of particles and is either still a leaf or empty
		if ((node->num_ > maxLeaf_) && (node->depth_ != maxDepth_)) {  // If node is no longer a leaf node
			node->leaf_ = false;  // node was previously a leaf node
			build(node, node->depth_);  // build child nodes from node
		}
		else {  // If node is a leaf node or an empty node
			if (node->num_ == 0) {  // empty case
				node->x_ = 0; node->y_ = 0; node->z_ = 0;
				delete[] node->Objs_; node->Objs_ = nullptr;
				node->leaf_ = false;  // node was previously a leaf node
			}
			else {  // leaf case
				vec3* V = new vec3[node->num_];  // get c.o.d. coords
				node->x_ = 0; node->y_ = 0; node->z_ = 0;
				for (int j = 0; j < node->num_; j++) {  // Calculate c.o.d. of node by producing array of vec3=length of object array
					V[j] = func(node->Objs_[j]);
					node->x_ += V[j].x_;
					node->y_ += V[j].y_;
					node->z_ += V[j].z_;
				}
				delete[] V;
				// set c.o.d.
				node->setX(node->x_ / node->num_);
				node->setY(node->y_ / node->num_);
				node->setZ(node->z_ / node->num_);
			}
		}
	}
	else if (node->num_ > 0) {  // If node previously was a parent node
								// Node has lost enough particles to become a leaf node or empty
								// Node's children have updated and node needs updating
		int oldNum = node->num_;
		int newNum = 0;  // Calculate new node->num_
		for (int i = 0; i < 8; i++) {
			if (node->child_[i]) {
				newNum += node->child_[i]->num_;
			}
		}
		if (oldNum == newNum) {  // The node and its children are all updated
			return;
		}
		// Node needs updating
		node->num_ = newNum;
		if ((node->num_ > maxLeaf_) && (node->depth_ != maxDepth_)) {  // If node is a parent node
			for (int i = 0; i < 8; i++) {  // Calculate c.o.d. of Node by iterating over new children statistics
				node->x_ += (node->child_[i]->x_ * node->child_[i]->num_);
				node->y_ += (node->child_[i]->y_ * node->child_[i]->num_);
				node->z_ += (node->child_[i]->z_ * node->child_[i]->num_);
			}
			node->setX(node->x_ / node->num_); node->setY(node->y_ / node->num_); node->setZ(node->z_ / node->num_);
		}
		else {  // If node has changed to a leaf node or an empty node
			if (node->num_ == 0) {  // empty case
				node->x_ = 0; node->y_ = 0; node->z_ = 0;
				node->Objs_ = nullptr;
				for (int i = 0; i < 8; i++) {
					delete node->child_[i];
					node->child_[i] = nullptr;
				}
			}
			else {  // leaf case
				Tname* objs = new Tname[node->num_];
				int counter = 0;
				for (int i = 0; i < 8; i++) {  // Iterate over nodes, passing objects to parent node
					for (int j = 0; j < node->child_[i]->num_; j++) {
						objs[counter] = node->child_[i]->Objs_[j];
						counter++;
					}
					delete node->child_[i];
					node->child_[i] = nullptr;
				}
				node->Objs_ = objs;
				vec3* V = new vec3[node->num_];  // get c.o.d. coords
				std::pair<vec3, double> result;
				node->x_ = 0; node->y_ = 0; node->z_ = 0;
				for (int j = 0; j < node->num_; j++) {  // Calculate c.o.d. of node by producing array of vec3=length of object array
					V[j] = func(node->Objs_[j]);
					node->x_ += V[j].x_;
					node->y_ += V[j].y_;
					node->z_ += V[j].z_;
				}
				node->leaf_ = true;  // Node is now a leaf node
				delete[] V;
				// set c.o.d.
				node->setX(node->x_ / node->num_);
				node->setY(node->y_ / node->num_);
				node->setZ(node->z_ / node->num_);
			}
		}
	}
	if (node->parent_) {  // Not at root
		updateNode(node->parent_);  // Repeat for parent node
	}
	return;
}

template <typename Tname>
template <moveable>
void Octree<Tname>::updateNode(Node<Tname>* node) {
	// Updates node statistics with information from object list
	if (node->leaf_) {  // If node was previously a leaf node
						// Node has gained enough particles and is no longer a leaf
						// Node has different number of particles and is either still a leaf or empty
		if ((node->num_ > maxLeaf_) && (node->depth_ != maxDepth_)) {  // If node is no longer a leaf node
			node->leaf_ = false;  // node was previously a leaf node
			build(node, node->depth_);  // build child nodes from node
		}
		else {  // If node is a leaf node or an empty node
			if (node->num_ == 0) {  // empty case
				node->x_ = 0; node->y_ = 0; node->z_ = 0;
				delete[] node->Objs_; node->Objs_ = nullptr;
				node->leaf_ = false;  // node was previously a leaf node
			}
			else {  // leaf case
				vec3* V = new vec3[node->num_];  // get c.o.d. coords
				node->x_ = 0; node->y_ = 0; node->z_ = 0;
				for (int j = 0; j < node->num_; j++) {  // Calculate c.o.d. of node by producing array of vec3=length of object array
					V[j] = func(node->Objs_[j]);
					node->x_ += V[j].x_;
					node->y_ += V[j].y_;
					node->z_ += V[j].z_;
				}
				delete[] V;
				// set c.o.d.
				node->setX(node->x_ / node->num_);
				node->setY(node->y_ / node->num_);
				node->setZ(node->z_ / node->num_);
			}
		}
	}
	else if (node->num_ > 0) {  // If node previously was a parent node
								// Node has lost enough particles to become a leaf node or empty
								// Node's children have updated and node needs updating
		int oldNum = node->num_;
		int newNum = 0;  // Calculate new node->num_
		for (int i = 0; i < 8; i++) {
			if (node->child_[i]) {
				newNum += node->child_[i]->num_;
			}
		}
		if (oldNum == newNum) {  // The node and its children are all updated
			return;
		}
		// Node needs updating
		node->num_ = newNum;
		if ((node->num_ > maxLeaf_) && (node->depth_ != maxDepth_)) {  // If node is a parent node
			for (int i = 0; i < 8; i++) {  // Calculate c.o.d. of Node by iterating over new children statistics
				node->x_ += (node->child_[i]->x_ * node->child_[i]->num_);
				node->y_ += (node->child_[i]->y_ * node->child_[i]->num_);
				node->z_ += (node->child_[i]->z_ * node->child_[i]->num_);
			}
			node->setX(node->x_ / node->num_); node->setY(node->y_ / node->num_); node->setZ(node->z_ / node->num_);
		}
		else {  // If node has changed to a leaf node or an empty node
			if (node->num_ == 0) {  // empty case
				node->x_ = 0; node->y_ = 0; node->z_ = 0;
				node->Objs_ = nullptr;
				for (int i = 0; i < 8; i++) {
					delete node->child_[i];
					node->child_[i] = nullptr;
				}
			}
			else {  // leaf case
				Tname* objs = new Tname[node->num_];
				int counter = 0;
				for (int i = 0; i < 8; i++) {  // Iterate over nodes, passing objects to parent node
					for (int j = 0; j < node->child_[i]->num_; j++) {
						objs[counter] = std::move(node->child_[i]->Objs_[j]);
						counter++;
					}
					delete node->child_[i];
					node->child_[i] = nullptr;
				}
				node->Objs_ = objs;
				vec3* V = new vec3[node->num_];  // get c.o.d. coords
				std::pair<vec3, double> result;
				node->x_ = 0; node->y_ = 0; node->z_ = 0;
				for (int j = 0; j < node->num_; j++) {  // Calculate c.o.d. of node by producing array of vec3=length of object array
					V[j] = func(node->Objs_[j]);
					node->x_ += V[j].x_;
					node->y_ += V[j].y_;
					node->z_ += V[j].z_;
				}
				node->leaf_ = true;  // Node is now a leaf node
				delete[] V;
				// set c.o.d.
				node->setX(node->x_ / node->num_);
				node->setY(node->y_ / node->num_);
				node->setZ(node->z_ / node->num_);
			}
		}
	}
	if (node->parent_) {  // Not at root
		updateNode(node->parent_);  // Repeat for parent node
	}
	return;
}

// Octree move/copy octree object functions
template <typename Tname>
template <moveable>
Tname* Octree<Tname>::moveTreeData(Node<Tname>* node, Tname* ObjArr, bool homeNode) {
	static int counter;
	if (homeNode) {
		homeNode = false;
		counter = 0;  // initialise array element counter
		int size = node->num_;
		if (ObjArr == nullptr) { ObjArr = new Tname[size](); }  // Create new array of size of Octree
		if (!!node->child_[0] + !!node->child_[1] + !!node->child_[2] + !!node->child_[3] +
			!!node->child_[4] + !!node->child_[5] + !!node->child_[6] + !!node->child_[7] == 0) {  // If no child nodes beyond homeNode
			for (int j = 0; j < node->num_; j++) {
				ObjArr[counter] = std::move(node->Objs_[j]);  // Move element to array
				counter++;  // increase iterator
			}
			node->num_ = 0;  // No objects in node->child_[i] now
			updateNode(node);  // Update node->child_[i] and its parent nodes
			return ObjArr;
		}
	}
	for (int i = 0; i < 8; i++) {  // Iterating over the child_ nodes
		if ((node->child_[i]->num_ > maxLeaf_) && (node->child_[i]->depth_ != maxDepth_)) {
			ObjArr = moveTreeData(node->child_[i], ObjArr, homeNode);  // repeat process with child_ nodes
																	  // returns array of children's particles
		}
		else {
			if (node->child_[i]->Objs_ != nullptr) {
				for (int j = 0; j < node->child_[i]->num_; j++) {
					ObjArr[counter] = std::move(node->child_[i]->Objs_[j]);  // Move element to array
					counter++;  // increase iterator
				}
				node->child_[i]->num_ = 0;  // No objects in node->child_[i] now
				updateNode(node->child_[i]);  // Update node->child_[i] and its parent nodes
			}
		}
	}
	return ObjArr;  // return array of objects
}

template <typename Tname>
template <copyable>
Tname* Octree<Tname>::copyTreeData(Node<Tname>* node, Tname* ObjArr, bool homeNode) const {
	static int counter;
	if (homeNode) {
		homeNode = false;
		counter = 0;  // initialise array element counter
		int size = node->num_;
		if (ObjArr == nullptr) { ObjArr = new Tname[size](); }  // Create new array of size of node->num_
		if (!!node->child_[0] + !!node->child_[1] + !!node->child_[2] + !!node->child_[3] +
			!!node->child_[4] + !!node->child_[5] + !!node->child_[6] + !!node->child_[7] == 0) {  // If no child nodes beyond homeNode
			for (int j = 0; j < node->num_; j++) {
				ObjArr[counter] = node->Objs_[j];  // Move element to array
				counter++;  // increase iterator
			}
			return ObjArr;
		}
	}
	for (int i = 0; i < 8; i++) {  // Iterating over the child_ nodes
		if ((node->child_[i]->num_ > maxLeaf_) && (node->child_[i]->depth_ != maxDepth_)) {
			ObjArr = copyTreeData(node->child_[i], ObjArr, homeNode);  // repeat process with child_ nodes
																	  // returns array of children's particles
		}
		else {
			if (node->child_[i]->Objs_ != nullptr) {
				for (int j = 0; j < node->child_[i]->num_; j++) {
					ObjArr[counter] = node->child_[i]->Objs_[j];  // Copy over element to array
					counter++;  // increase iterator
				}
			}
		}
	}
	return ObjArr;  // return array of objects
}

// Octree add to octree functions
template <typename Tname>
template<copyableOnly>
void Octree<Tname>::addToTree(Tname Obj) {
	vec3 coords = func(Obj);  // using function provided to get obj coords
	Node<Tname>* destNode = findLeafNode(coords.x_, coords.y_, coords.z_);
	if (!destNode) {  // object is out of bounds. Call addToTree(Tname*).
		Tname* ObjArr = new Tname[1]; ObjArr[0] = std::move(Obj);
		addToTree(ObjArr, 1);  // rebuild tree
	}
	else {  // object is in bounds. Update destination node and parents
		Tname* destObjs = new Tname[destNode->num_ + 1];
		for (int k = 0; k < destNode->num_; k++) {
			destObjs[k] = destNode->Objs_[k];
		}
		// Move object to leaf node, adjust array sizes appropriately
		destObjs[destNode->num_] = Obj;
		destNode->num_++;
		if (destNode->Objs_ != nullptr) {
			delete[] destNode->Objs_;
		}
		destNode->Objs_ = destObjs;
		updateNode(destNode);  // Update destination node, building more leaves if needed, update extended family nodes
	}
	return;
}

template <typename Tname>
template<moveable>
void Octree<Tname>::addToTree(Tname Obj) {
	vec3 coords = func(Obj);  // using function provided to get obj coords
	Node<Tname>* destNode = findLeafNode(coords.x_, coords.y_, coords.z_);
	if (!destNode) {  // object is out of bounds. Call addToTree(Tname*).
		Tname* ObjArr = new Tname[1]; ObjArr[0] = std::move(Obj);
		addToTree(ObjArr, 1);  // rebuild tree
	}
	else {  // object is in bounds. Update destination node and parents
		Tname* destObjs = new Tname[destNode->num_ + 1];
		for (int k = 0; k < destNode->num_; k++) {
			destObjs[k] = std::move(destNode->Objs_[k]);
		}
		// Move object to leaf node, adjust array sizes appropriately
		destObjs[destNode->num_] = std::move(Obj);
		destNode->num_++;
		if (destNode->Objs_ != nullptr) {
			delete[] destNode->Objs_;
		}
		destNode->Objs_ = destObjs;
		updateNode(destNode);  // Update destination node, building more leaves if needed, update extended family nodes
	}
	return;
}

template <typename Tname>
template<copyableOnly>
void Octree<Tname>::addToTree(Tname* ObjArr, int size) {
	int treeTot = root_->num_ + size;
	Tname* Objects = copyTreeData(root_);  // get data
	root_->num_ = treeTot;  // add to root size
	Tname* allObjects = new Tname[treeTot];
	for (int j = 0; j < root_->num_ - size; j++) {
		allObjects[j] = Objects[j];
	}
	for (int j = 0; j < size; j++) {
		allObjects[root_->num_ - size + j] = ObjArr[j];  // copy over array of objects
	}
	for (int i = 7; i >= 0; i--) {
		if (root_->child_[i]) { delete root_->child_[i]; }  // delete children
		root_->child_[i] = nullptr;
	}
	delete[] Objects;
	root_->Objs_ = allObjects;
	root_->x_ = 0; root_->y_ = 0; root_->z_ = 0;  // reinitialise root C.o.D.
	build(root_);  // rebuild tree
	return;
}

template <typename Tname>
template<moveable>
void Octree<Tname>::addToTree(Tname* ObjArr, int size) {
	int treeTot = root_->num_ + size;
	Tname* Objects = moveTreeData(root_);  // get data
	root_->num_ = treeTot;  // add to root size
	Tname* allObjects = new Tname[treeTot];
	for (int j = 0; j < root_->num_ - size; j++) {
		allObjects[j] = std::move(Objects[j]);
	}
	for (int j = 0; j < size; j++) {
		allObjects[root_->num_ - size + j] = std::move(ObjArr[j]);  // copy over array of objects
	}
	for (int i = 7; i >= 0; i--) {
		if (root_->child_[i]) { delete root_->child_[i]; }  // delete children
		root_->child_[i] = nullptr;
	}
	delete[] Objects;
	root_->Objs_ = std::move(allObjects);
	root_->x_ = 0; root_->y_ = 0; root_->z_ = 0;  // reinitialise root C.o.D.
	build(root_);  // rebuild tree
	return;
}

#endif // !OCTREE_H