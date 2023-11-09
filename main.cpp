#include <iostream>
#include <vector>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <iomanip>


/*
Calculate the square of the distance between two points.
*/
template <typename T>
T distance_squared(const T* point_a, const T* point_b, size_t N) {
	T d_squared = T();
	for (size_t i = 0; i < N; ++i) {
		auto delta = point_a[i] - point_b[i];
		d_squared += delta * delta;
	}
	return d_squared;
}

template <typename T>
T distance(const T* point_a, const T* point_b, size_t N) {
	return std::sqrt(distance_squared(point_a, point_b,N));
}

/*
Here is your parameters of E2 LSH:
L (the number of hash funcitons) = 233.913
T (collision threshold) = 113.713
R (query range) = 0.565685
Alpha (collision percentage) = 0.486135
W (bucket width) [user input] =  1
C (approximation ratio) [user input] =  1.5
Delta (false negative rate) [user input] =  0.05
Theta (false negative rate) [user input] =  0.3
K (# projections in each hash function) [user input] =  1
D Dimentions
*/

unsigned hashf(float x)
{
    union
    {
        float f;
        unsigned u;
    };
    f = x;
    u ^= u >> 13;
    u *= 0x5bd1e995;
    return u ^ (u >> 15);    
}


class LSH {
    public:
    using vector3Df = std::vector<std::vector<std::vector<float>>>;
    using vector1Df = std::vector<float>;
    LSH() {
        LSH_init();
    }
    LSH( int L, int K, int D ) : L(L), K(K), D(D) {
        LSH_init();
    }

    double calc_distance(const float * first, const float *second) {
        return ::distance(first,second,D);
    }

    void print(const float * vec) {
        if (vec == nullptr) {
            std::cout << "null";
            
        }
        else
        for (int i = 0; i < D; ++i)
        {
            std::cout << std::fixed << std::setprecision(4) << std::setw(8) << vec[i];            
        }
        std::cout << std::endl;
    }

    std::unordered_set<unsigned> computeHash( const float * vec)
    {
            //loop through # of hash function group
            std::unordered_set<unsigned> vL;

            for (int l = 0; l < L; ++l) {
                float hashFinalValue = 0;
                //loop through the inner of a hash function group
                for (int k = 0; k < K; ++k) {
                    float dTemp = 0;
                    //loop through all the dimensions
                    for (int d = 0; d < D; ++d) {
                        //vector(math) multiply to make projection
                        dTemp += (randomLine)[l][k][d]* vec[d];
                    }
                    //assign hashvalue into bucket
                    float hashvalue =  floor(dTemp/1);
                    //std::cout << int(hashvalue) << "|";
                    //merge hash group results **see documentation
                    hashFinalValue = hashvalue*randomVector[k] + hashFinalValue;
                }
                //vL.push_back(hashFinalValue);
                vL.insert( hashf( hashFinalValue ) );
            }
            return vL;
    }

    void insert( const float * vec )
    {
        for (auto hash : computeHash(vec)) {
            hashed.insert( std::make_pair( hash, vec));
            //std::cout << hash << " ";
        }
    }

    const float * closest( const float * vec) 
    {
        const float * result = nullptr;
        for (auto hash : computeHash(vec) ) 
        {             
            auto range = hashed.equal_range(hash);
            
            float distance = 9999999;
            for (auto it = range.first; it != range.second; ++it)
            {
               auto dis = distance_squared( vec, it->second, D );
               if (dis < distance) {
                   distance = dis;
                   result = it->second;
               }
            }
        }
        return result;        
    }
    
    const float * closest( const float * vec, float break_distance) 
    {
        const float * result = nullptr;
        break_distance *= break_distance;
        for (auto hash : computeHash(vec) ) 
        {             
            auto range = hashed.equal_range(hash);
            
            float distance = 9999999;
            for (auto it = range.first; it != range.second; ++it)
            {
               auto dis = distance_squared( vec, it->second, D );
               if (dis < break_distance) {
                   return it->second;
               }
            }
        }
        return nullptr;        
    }


    const auto& getHashed() {
        return hashed;
    }

    private:

    void LSH_init() {
        generateRandomLine();
        generateUniformRandomVector();

    }


    /**
     * This function generates the random projected lines, it uses normal distribution
     * @return the 3D vector of randomLine generated
     */
    void generateRandomLine(){

        std::random_device rd;
        std::mt19937 gen(rd());

        for (int i = 0; i < L; i++) {
            auto & vK  = randomLine.emplace_back();
            for (int j = 0; j < K; ++j) {
                //generate random according normal distribution
                std::normal_distribution<float> distribution(0.5,0.5);
                auto & vD = vK.emplace_back();
                for (int k = 0; k < D; ++k) {
                    auto d = distribution(gen);
                    //std::cout << "{" << d << "} ";
                    vD.push_back(d);
                }
            }
        }
    }

    /**
     * This funciton generates k number of random vectors, it uses uniform random generator
     * @param number number of randomNumber, usually = K
     * @param maxium the max value for the random generation
     * @return 1D vector of random numbers
     */
    void generateUniformRandomVector(){
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0, 1);

        for (int i = 0; i < K; ++i) {
            randomVector.push_back(dis(gen));
        }

    }

    int L=8,K=1,D=16;
    vector3Df randomLine;
    vector1Df randomVector;
    std::unordered_multimap<float,const float*> hashed;
};


int main() {

    
    std::vector<std::vector<float>> dataset(32, std::vector<float>(16));
    LSH lsh(8,1,dataset[0].size());
    float x = 0;
    float y = 0.001;

    for (auto & row : dataset)
    {      
        
        for (auto & v : row)  {
            v = sin(x) * y;
            x += 6.28f / row.size();            
        }

        lsh.print(&row[0]);
        
        auto v = lsh.closest(&row[0], 0.009);

        lsh.print(v);

        if ( v != nullptr ) 
        {
           std::cout << "distance=" << lsh.calc_distance(&row[0], v);
        } else {
            lsh.insert(&row[0]);
        }

        std::cout << std::endl;
    
        
        y += y;
        if (y > 1)
            y = 1;

    }




}
