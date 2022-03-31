// Code tests
// Compile with g++ -std=c++20 -o EllieOctreeTests -pthread -O3 tests.cpp
// Run with ./EllieOctreeTests

#include "../include/EllieOctree.hpp"
#include <chrono>
#include <random>
#include <memory>

vec3 vec3SmrtPntrCoords(std::shared_ptr<vec3>& vec) {
	vec3 V;
	V.x_ = vec.get()->x_; V.y_ = vec.get()->y_; V.z_ = vec.get()->z_;
	return V;
}  // returns vec3 obj coords (for octree)

vec3 vec3Coords(vec3& vec) {
	vec3 V;
	V.x_ = vec.x_; V.y_ = vec.y_; V.z_ = vec.z_;
	return V;
}  // returns vec3 obj coords (for octree)


int main() {
	int numParts = 1000000;
	// Initialise random generator
	std::random_device rd;
	std::mt19937 gen(rd());  // create random numbers
	std::uniform_real_distribution<> uni(0, 1);
	// Make list of physical particles
	vec3* particleList = new vec3[numParts];
	double maxX1 = 0;
	double maxY1 = 0;
	double maxZ1 = 0;
	for (int i = 0; i < numParts; i++) {
		double x = uni(gen); maxX1 = (x > maxX1) ? x : maxX1;
		double y = uni(gen); maxY1 = (y > maxY1) ? y : maxY1;
		double z = uni(gen); maxZ1 = (z > maxZ1) ? z : maxZ1;
		particleList[i] = vec3{ x, y, z };
	}
	// Make octree of physical particles
	Octree<vec3>* octree3 = new Octree<vec3>(particleList, vec3Coords, numParts, 3, 16, 0, 1, 0, 1, 0, 1);
	// Get a node, regenerate particle positions
	Node<vec3>* leafNode = octree3->findLeafNode(maxX1 - 0.0001, maxY1 - 0.0001, maxZ1 - 0.0001);
	for (int i = 0; i < leafNode->num_; i++) {
		leafNode->Objs_[i].x_ = uni(gen);
		leafNode->Objs_[i].y_ = uni(gen);
		leafNode->Objs_[i].z_ = uni(gen);
	}
	// Update octree
	octree3->updateTree(octree3->getRoot());
	// Memory cleanup
	delete octree3;
	for (int i = 0; i < numParts; i++) {
		particleList[i].~vec3();
	}
	delete[] particleList;

	// Make list of pointers to particles
	//Particle** particles = new Particle*[numParts]();
	std::shared_ptr<vec3>* particles0 = new std::shared_ptr<vec3> [numParts]();
	double maxX = 0;
	double maxY = 0;
	double maxZ = 0;
	for (int i = 0; i < numParts; i++) {
		double x = uni(gen); maxX = (x > maxX) ? x : maxX;
		double y = uni(gen); maxY = (y > maxY) ? y : maxY;
		double z = uni(gen); maxZ = (z > maxZ) ? z : maxZ;
		std::shared_ptr<vec3> particle(new vec3{ x, y, z });
		particles0[i] = std::move(particle);
	}
	std::cout << numParts << " particles constructed." << std::endl;
	// Build an octree
	auto t1a = std::chrono::high_resolution_clock::now();
	Octree<std::shared_ptr<vec3>>* octree = new Octree<std::shared_ptr<vec3>>(particles0, vec3SmrtPntrCoords, numParts, 3, 16, 0, 1, 0, 1, 0, 1);		// 4ms for 10k particles
	auto t2a = std::chrono::high_resolution_clock::now();
	std::cout << "Octree build time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2a - t1a).count() << "ms" << std::endl;

	std::shared_ptr<vec3>* particles = octree->copyTreeData(octree->getRoot());
	//delete[] parts;

	// Change 1 particle position slightly, and update octree
	// Edit particle
	std::cout << particles[0].get()->x_ << std::endl;
	particles[0].get()->x_ = particles[0].get()->x_ + 0.05;
	std::cout << particles[0].get()->x_ << std::endl;
	// Test time to retrieve octree data
	auto t1b = std::chrono::high_resolution_clock::now();
	Node<std::shared_ptr<vec3>>* root = octree->getRoot();
	std::cout << "Retrieved root." << std::endl;
	std::shared_ptr<vec3>* octreeParticles = octree->copyTreeData(octree->getRoot());  // collect octree particles
	auto t2b = std::chrono::high_resolution_clock::now();
	std::cout << "Octree data retrieval time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2b - t1b).count() << "ms" << std::endl;  // 10us for 10k particles
	//octree->addToTree(octreeParticles, numParts);
	delete[] octreeParticles;

	// Update the octree
	auto t3a = std::chrono::high_resolution_clock::now();
	std::shared_ptr<vec3>* lostObjs1 = octree->updateTree(octree->getRoot());
	auto t4a = std::chrono::high_resolution_clock::now();
	std::cout << "Update octree time (1 changed particle): " << std::chrono::duration_cast<std::chrono::milliseconds>(t4a - t3a).count() << "ms" << std::endl;
	if (lostObjs1) { delete[] lostObjs1; }

	// Change particle positions slightly, and update octree
	// Edit particles
	for (int i = 0; i < numParts; i++) {
		particles[i].get()->x_ += 0.1;
	}
	// Update the octree
	auto t3b = std::chrono::high_resolution_clock::now();
	std::shared_ptr<vec3>* lostObjs2 = octree->updateTree(octree->getRoot());
	auto t4b = std::chrono::high_resolution_clock::now();
	int lostParticles = numParts - octree->getRoot()->num_;
	std::cout << "Update octree time (" << numParts << " changed particles): " << std::chrono::duration_cast<std::chrono::milliseconds>(t4b - t3b).count() << "ms" << std::endl;
	std::cout << numParts << ", " << octree->getRoot()->num_ << std::endl;
	std::cout << "Lost " << lostParticles << " particles." << std::endl; // We should lose ~10% of the particles

	// Reset lost particles to old positions
	for (int i = 0; i < lostParticles; i++) {
		lostObjs2[i].get()->x_ = lostObjs2[i].get()->x_ - 0.1;
	}
	// Add lost particles to the tree
	auto t5 = std::chrono::high_resolution_clock::now();
	if (lostObjs2) { octree->addToTree(lostObjs2, lostParticles); }
	auto t6 = std::chrono::high_resolution_clock::now();
	std::cout << "Time to add " << lostParticles << " particles: " << std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count() << "ms" << std::endl;
	delete[] lostObjs2;
	// Retrieve particles, and check the tree has updated correctly
	std::cout << "Octree root node now has " << octree->getRoot()->num_ << " particles." << std::endl;
	std::shared_ptr<vec3>* newOctreeParticles = octree->copyTreeData(octree->getRoot());
	std::cout << "Particle number " << numParts << " has x-coordinate: " << newOctreeParticles[numParts - 1]->x_ << std::endl;
	delete[] newOctreeParticles;

	// Regenerate particle positions, and update octree
	// Edit particles with massively changed positions
	for (int i = 0; i < numParts; i++) {
		particles[i].get()->x_ = uni(gen);
		particles[i].get()->y_ = uni(gen);
		particles[i].get()->z_ = uni(gen);
	}
	auto t7 = std::chrono::high_resolution_clock::now();
	octree->updateTree(octree->getRoot());
	auto t8 = std::chrono::high_resolution_clock::now();
	std::cout << "Octree update time with massive particle changes: " << std::chrono::duration_cast<std::chrono::milliseconds>(t8 - t7).count() << "ms" << std::endl;

	// Destroy and rebuild octree entirely
	auto t9 = std::chrono::high_resolution_clock::now();
	delete octree;
	Octree<std::shared_ptr<vec3>>* octree2 = new Octree<std::shared_ptr<vec3>>(particles, vec3SmrtPntrCoords, numParts, 3, 16, 0, 1, 0, 1, 0, 1);
	auto t10 = std::chrono::high_resolution_clock::now();
	std::cout << "Time to destroy and rebuild octree (with data already retrieved): " << std::chrono::duration_cast<std::chrono::milliseconds>(t10 - t9).count() << "ms" << std::endl;
	std::cout << octree2->getRoot()->xMin_ << ", " << octree2->getRoot()->xMax_ << std::endl;

	// Make a new particle in the box limits and add it
	double newX1 = uni(gen);
	double newY1 = uni(gen);
	double newZ1 = uni(gen);
	std::shared_ptr<vec3> newParticle1(new vec3{ newX1, newY1, newZ1 });
	auto t11 = std::chrono::high_resolution_clock::now();
	octree2->addToTree(std::move(newParticle1));
	auto t12 = std::chrono::high_resolution_clock::now();
	std::cout << "Time to add particle inside octree bounds (in microseconds): " << std::chrono::duration_cast<std::chrono::microseconds>(t12 - t11).count() << "us" << std::endl;

	// Make a new particle outside the box limits and add it (should be == destroy+rebuild)
	std::shared_ptr<vec3> newParticle2(new vec3{ 1.1, 1.1, 1.1 });
	auto t13 = std::chrono::high_resolution_clock::now();
	octree2->addToTree(std::move(newParticle2));
	auto t14 = std::chrono::high_resolution_clock::now();
	std::cout << "Time to add particle outside of octree bounds (in milliseconds): " << std::chrono::duration_cast<std::chrono::milliseconds>(t14 - t13).count() << "ms" << std::endl;

	std::cout << octree2->getDataSize(octree2->getRoot()) << std::endl;
	std::cout << octree2->getDataSize(octree2->getRoot()->child_[0]) << std::endl;
	std::cout << octree2->getRoot()->xMin_ << ", " << octree2->getRoot()->xMax_ << std::endl;

	delete octree2;
	newParticle1.reset();
	newParticle2.reset();

	return 0;
}