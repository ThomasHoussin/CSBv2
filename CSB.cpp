#pragma GCC optimize "O3,omit-frame-pointer,inline,fast-math"

//collision algorithm based on http://files.magusgeek.com/csb/csb_en.html

#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <random>
#include <memory>
#include <chrono>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//default max time allowed for IA computing (in ms)
#define DEFAULT_MAX_TIME 145
//number of moves per pod in a solution
#define LENGTH 6
//max distance squared /100 on map, used for normalization
#define MAX_DIST 3610000.0
//shield activation probability
#define ASHIELD 0.1
//probability for max and min thrust
#define MAXTHRUSTP 0.15
#define MINTHRUSTP 0.15
//constant for collisions
#define EPSILON 0.00001

//exponential function approximation
inline
double exp(double x) {
	  x = 1.0 + x / 16.0;
	  x *= x; x *= x; x *= x; x *= x;
	  return x;
}

//singleton class to handle random numbers
class Random {
private:
    std::mt19937 engine ;
    std::uniform_real_distribution<double> fdistribution ;
    std::uniform_int_distribution<int> angle_distribution ;
    std::uniform_int_distribution<int> thrust_distribution ;
    std::uniform_int_distribution<int> length_distribution ;

    //copies are not allowed
    Random(Random const&) {}
	void operator=(Random const&) {}

	//constructor is private
	//getInstance is the way to get a Random object (unique)
	Random() {
	    static std::random_device rd; // obtain a random number from hardware
	    static std::seed_seq seed{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()}; //build seed sequence
	    engine = std::mt19937(seed); // seed the generator

	    fdistribution = std::uniform_real_distribution<double>(0.0,1.0) ;
	    angle_distribution = std::uniform_int_distribution<int>(-9,9) ;
	    thrust_distribution = std::uniform_int_distribution<int>(-20, 20) ;
	    length_distribution = std::uniform_int_distribution<int>(0, 2 * LENGTH - 1) ;
	}

public:
	~Random(){}

	static Random& getInstance() {
		static Random instance ;
		return instance ;
	}
	double nextFloat() {
		return fdistribution(engine) ;
	}
	int nextRAngle() {
		return angle_distribution(engine);
	}
	int nextRThrust() {
		//use of 5 * random(0,40) should reduce the search space without much effect on the result
		return thrust_distribution(engine) * 5;
	}
	int nextRLength() {
		return length_distribution(engine);
	}
	unsigned long nextRandom() {
		return engine() ;
	}
};

class Move {
public:
	Move(double angle = 0., int thrust=100) :
		angle(angle), thrust(thrust)
	{ }

	double getAngle() const {
		return angle;
	}

	void setAngle(double angle) {
		this->angle = angle;
	}

	int getThrust() const {
		return thrust;
	}

	void setThrust(int thrust) {
		this->thrust = thrust;
	}

	//generate a neighbor move : the basic idea is to have a better neighborhood than just a random move, to help SA convergence
	//we use a fixed probability for shield activation for the turn
	//for thrust  : fixed probability for min and max thrust (fast chage is required for exemple when crossing CP)
	//in other case, we use a random value ; this value is added to the current thrust
	//for example : if thrust is 100, random value between -100 and 100, result will be between 0 and 200
	//if thrust is 200, result will be between 100 and 200, with increased probability on 200
	void mutate() {
		Random &r = Random::getInstance() ;
		double p = r.nextFloat() ;
		if(p < ASHIELD) {
			this->thrust = -1 ;
		}
		else if(p < ASHIELD + MINTHRUSTP) {
		    thrust = 0 ;
		}
		else if(p < ASHIELD + MINTHRUSTP + MAXTHRUSTP) {
		    thrust = 200 ;
		}
		else {
			thrust += r.nextRThrust() ;
			if(thrust > 200) thrust = 200 ;
			else if(thrust < 0) thrust = 0 ;
		}

		angle += r.nextRAngle() ;
		if(angle > 18) angle = 18 ;
		else if(angle < -18) angle = -18 ;
	}

	static Move randomMove() {
		Move m ;
		m.mutate() ;
		return m ;
	}

private :
	double angle; // Between -18.0 and +18.0
    int thrust; // Between -1 and 200
};

std::ostream& operator <<(std::ostream& Stream, const Move& move)
{
    Stream << "Angle : " << move.getAngle() << " Thrust : " << move.getThrust() ;
    return Stream; // N'oubliez pas de renvoyer le flux, afin de pouvoir chaîner les appels
}

class Solution {
private:
	double score ;
	int savedPosition ;
	double savedScore ;
	Move saved ;

public:
	//LENGTH first elements : first pod
	//LENGTH last elements : second pod
	Move solution[2 * LENGTH] ;
	Solution& operator=(Solution const& s) ;

	void save(int n) {
		savedPosition = n ;
		savedScore = score ;
		saved = solution[n] ;
	}

	void restore() {
		solution[savedPosition] = saved ;
		score = savedScore ;
		savedPosition = -1 ;
	}

	double getSavedScore() const {
		return savedScore ;
	}

	int getSavedPosition() const {
		return savedPosition ;
	}

	Solution(bool buildRandom = false) : score(std::numeric_limits<double>::lowest()), savedPosition(-1), savedScore(0.), saved(Move(0,0)) {
		if(buildRandom) {
			for(int i = 0 ; i < 2* LENGTH ; ++i) {
				this->solution[i] = Move::randomMove() ;
			}
		}
		else {
			for(int i = 0 ; i < 2* LENGTH ; ++i) {
				this->solution[i] = Move() ;
			}
		}
	}

	Solution(const Solution& s) : score(s.score), savedPosition(s.savedPosition), savedScore(s.savedScore), saved(s.saved) {
		std::copy(std::begin(s.solution), std::end(s.solution), std::begin(this->solution));
	}

	void shiftLeft() {
		std::copy(std::begin(this->solution)+1, std::end(this->solution),std::begin(this->solution)) ;
		this->solution[LENGTH - 1] = Move() ;
	}

	static Solution* randomSolution() {
		Solution* s = new Solution() ;
		for(int i = 0 ; i < 2* LENGTH ; ++i) {
			s->solution[i] = Move::randomMove() ;
		}
		return s ;
	}

	void mutate(){
		Random& r = Random::getInstance() ;
		int l = r.nextRLength() ;
		this->save(l);
		solution[l].mutate() ;
	}

	double getScore() const {
		return score;
	}

	void setScore(double score) {
		this->score = score;
	}
};

Solution& Solution::operator=(Solution const& s){
	this->score = s.score ;
	std::copy(std::begin(s.solution), std::end(s.solution), std::begin(this->solution)) ;
	return *this ;
}

class Point {
public:
	double x ;
	double y ;

	Point() : x(0.), y(0.){
	}
	Point(double x, double y): x(x), y(y) {
	}
	double distance(double x, double y) const {
		return sqrt((this->x - x) * (this->x - x) + (this->y - y) * (this->y - y)) ;
	}
	double distance(const Point &p) {
		return this->distance(p.x, p.y) ;
	}
	double distanceSq(double x, double y) const {
		return (this->x - x) * (this->x - x) + (this->y - y) * (this->y - y) ;
	}
	double distanceSq(const Point &p) {
		return this->distanceSq(p.x, p.y) ;
	}

	double getX() const {
		return x;
	}

	void setX(double x) {
		this->x = x;
	}

	double getY() const {
		return y;
	}

	void setY(double y) {
		this->y = y;
	}

	//It allows us to find the closest point on a line (described by two points) to another point.
	Point closest(const Point &a, const Point &b) {
	    double da = b.y - a.y;
	    double db = a.x - b.x;
	    double c1 = da*a.x + db*a.y;
	    double c2 = -db*this->x + da*this->y;
	    double det = da*da + db*db;
	    double cx = 0;
	    double cy = 0;

	    if (det != 0) {
	        cx = (da*c1 - db*c2) / det;
	        cy = (da*c2 + db*c1) / det;
	    } else {
	        // The point is already on the line
	        cx = this->x;
	        cy = this->y;
	    }

	    return Point(cx, cy);
	}

protected:

};

//forward declarations
class Unit ;
class Pod ;

class Collision {
public:
	Unit* a ;
	Unit* b;
	double t ;

	Collision(Unit* a, Unit* b, double t) :
		a(a), b(b), t(t) {

	}

private:

};

class Unit:public Point {
public:
	Unit(double x, double y):
		Point(x, y),r(0), vx(0), vy(0)
	{ }

	virtual ~Unit() {}

	Unit(double x, double y, double vx, double vy, double r):
		Point(x, y), r(r), vx(vx), vy(vy)
	{ }

	double getR() const {
		return r;
	}

	double getVx() const {
		return vx;
	}

	void setVx(double vx) {
		this->vx = vx;
	}

	double getVy() const {
		return vy;
	}

	void setVy(double vy) {
		this->vy = vy;
	}

	Collision* collision(Unit* u) {
		// Square of the distance
	    double dist = this->distanceSq(*u);

	    // Sum of the radii squared
	    double sr = (this->r + u->r)*(this->r + u->r);

	    if (dist < sr) {
	        // Objects are already touching each other. We have an immediate collision.
	        return new Collision(this, u, 0.0);
	    }

	    // Optimisation. Objects with the same speed will never collide
	    if (this->vx == u->vx && this->vy == u->vy) {
	        return nullptr;
	    }

	    // We place ourselves in the reference frame of u. u is therefore stationary and is at (0,0)
	    double x = this->x - u->x;
	    double y = this->y - u->y;
	    Point myp(x, y);
	    double vx = this->vx - u->vx;
	    double vy = this->vy - u->vy;
	    Point up(0., 0.) ;

	    // We look for the closest point to u (which is in (0,0)) on the line described by our speed vector
	    Point p = up.closest(myp, Point(x + vx, y + vy));

	    // Square of the distance between u and the closest point to u on the line described by our speed vector
	    double pdist = up.distanceSq(p);

	    // Square of the distance between us and that point
	    double mypdist = myp.distanceSq(p);

	    // If the distance between u and this line is less than the sum of the radii, there might be a collision
	    if (pdist < sr) {
	     // Our speed on the line
	        double length = sqrt(vx*vx + vy*vy);

	        // We move along the line to find the point of impact
	        double backdist = sqrt(sr - pdist);
	        p.x = p.x - backdist * (vx / length);
	        p.y = p.y - backdist * (vy / length);

	        // If the point is now further away it means we are not going the right way, therefore the collision won't happen
	        if (myp.distanceSq(p) > mypdist) {
	            return nullptr;
	        }

	        pdist = p.distance(myp);

	        // The point of impact is further than what we can travel in one turn
	        if (pdist > length) {
	            return nullptr;
	        }

	        // Time needed to reach the impact point
	        double t = pdist / length;

	        return new Collision(this, u, t);
	    }

	    return nullptr;
	}


protected:
	double r ;
	double vx ;
	double vy ;
} ;



class Pod: public Unit {
public:
	Pod(double x = -1., double y=-1., double vx=0., double vy=0., double angle = -1, int nextCheckpointId = -1):
		Unit(x,y,vx,vy,400.0), angle(angle), nextCheckpointId(nextCheckpointId), timeout(100), shield(0),
		checked(0), hasBoost(true), finishTime(0)
	{ }

	double getAngle() const {
		return angle;
	}

	void setAngle(double angle) {
		this->angle = angle;
	}

	int getNextCheckpointId() const {
		return nextCheckpointId;
	}

	void setNextCheckpointId(int nextCheckpointId) {
		this->nextCheckpointId = nextCheckpointId;
	}

	bool isShieldActive() const {
		return this->shield == 4 ;
	}

	void activeShield() {
		this->shield = 4 ;
	}

	void decShield() {
		if(this->shield > 0) this->shield -- ;
	}

	void decTimeout() {
		this->timeout -- ;
	}

	int getTimeout() const {
		return this->timeout ;
	}

	void addFinishTime(double t) {
		this->finishTime += t ;
	}

	double getFinishTime() const {
		return this->finishTime ;
	}

	void updatePod(double x, double y, double vx, double vy, double angle, int nextCheckPointId) {
		this->x = x ;
		this->y = y ;
		this->vx = vx ;
		this->vy = vy ;
		this->angle = angle ;
		this->nextCheckpointId = nextCheckPointId ;
		this->finishTime = 0 ;
	}

	void checkedCP(int nNextCP) {
		this->timeout = 100 ;
		++ this->checked ;
		this->nextCheckpointId = nNextCP ;
	}

	//compute the angle that we should have to target a point
	double computeAngle(const Point& p) {
	    double d = this->distance(p);
	    double dx = (p.x - this->x) / d;
	    double dy = (p.y - this->y) / d;

	    // Simple trigonometry. We multiply by 180.0 / PI to convert radiants to degrees.
	    double a = acos(dx) * 180.0 / M_PI;

	    // If the point I want is below me, I have to shift the angle for it to be correct
	    if (dy < 0) {
	        a = 360.0 - a;
	    }
	    return a;
	}

	//return the difference between the target angle and the current angle
	double diffAngle(const Point& p) {
	    double a = this->computeAngle(p);

	    // To know whether we should turn clockwise or not we look at the two ways and keep the smallest
	    // The ternary operators replace the use of a modulo operator which would be slower
	    double right = this->angle <= a ? a - this->angle : 360.0 - this->angle + a;
	    double left = this->angle >= a ? this->angle - a : this->angle + 360.0 - a;

	    if (right < left) {
	        return right;
	    } else {
	        // We return a negative angle if we must rotate to left
	        return -left;
	    }
	}

	double diffSpeedAngle(const Point& p) {
		if(this->vx == 0 && this->vy == 0) return 0. ;
		//target angle
		double a = this->computeAngle(p) ;
		//angle given by speed vector
		double sa = this->computeAngle(Point(x + this->vx, y + this->vy)) ;

	    // To know whether we should turn clockwise or not we look at the two ways and keep the smallest
	    // The ternary operators replace the use of a modulo operator which would be slower
	    double right = sa <= a ? a - sa : 360.0 - sa + a;
	    double left = sa >= a ? sa - a : sa + 360.0 - a;

	    if (right < left) {
	        return right;
	    } else {
	        // We return a negative angle if we must rotate to left
	        return -left;
	    }
	}

	void rotate(const Point& p) {
	    double a = this->diffAngle(p);

	    // Can't turn by more than 18° in one turn
	    if (a > 18.0) {
	        a = 18.0;
	    } else if (a < -18.0) {
	        a = -18.0;
	    }

	    this->angle += a;

	    // The % operator is slow. If we can avoid it, it's better.
	    if (this->angle >= 360.0) {
	        this->angle = this->angle - 360.0;
	    } else if (this->angle < 0.0) {
	        this->angle += 360.0;
	    }
	}

	void rotate(double angle) {
		if(angle > 18.0) angle = 18.0 ;
		else if(angle < -18.0) angle = -18.0 ;

		this->angle += angle;

		// The % operator is slow. If we can avoid it, it's better.
		if (this->angle >= 360.0) {
		    this->angle = this->angle - 360.0;
		} else if (this->angle < 0.0) {
		    this->angle += 360.0;
		}
	}

	void boost(int thrust) {
	    if(thrust == -1) {
	    	this->activeShield() ;
	    	return;
	    }
	    else if(thrust == 650) {
	    	if(!hasBoost) thrust = 200 ;
	    	this->hasBoost = false ;
	    }
		//a pod which has activated its shield cannot accelerate for 3 turns
	    if (this->shield > 0) {
	        return;
	    }

	    // Conversion of the angle to radiants
	    double ra = this->angle * M_PI / 180.0;

	    // Trigonometry
	    this->vx += cos(ra) * thrust;
	    this->vy += sin(ra) * thrust;
	}

	void move(double t) {
	    this->x += this->vx * t;
	    this->y += this->vy * t;
	}

	void endTurn() {
	    this->x = round(this->x);
	    this->y = round(this->y);
	    this->vx = std::trunc(this->vx * 0.85);
	    this->vy = std::trunc(this->vy * 0.85);

	    // Don't forget that the timeout goes down by 1 each turn. It is reset to 100 when you pass a checkpoint
	    this->decTimeout() ;
	    this->decShield() ;
	}

	void bounce(Pod* u) {
		// If a pod has its shield active its mass is 10 otherwise it's 1
		double m1 = this->isShieldActive() ? 10 : 1;
		double m2 = u->isShieldActive() ? 10 : 1;
		double mcoeff = (m1 + m2) / (m1 * m2);

		double nx = this->x - u->x;
		double ny = this->y - u->y;

		// Square of the distance between the 2 pods.
		double nxnysquare = nx*nx + ny*ny;

		double dvx = this->vx - u->vx;
		double dvy = this->vy - u->vy;

		// fx and fy are the components of the impact vector. product is just there for optimisation purposes
		double product = nx*dvx + ny*dvy;
		double fx = (nx * product) / (nxnysquare * mcoeff);
		double fy = (ny * product) / (nxnysquare * mcoeff);

		// We apply the impact vector once
		this->vx -= fx / m1;
		this->vy -= fy / m1;
		u->vx += fx / m2;
		u->vy += fy / m2;

		// If the norm of the impact vector is less than 120, we normalize it to 120
		double impulse = sqrt(fx*fx + fy*fy);
		if (impulse < 120.0) {
			fx = fx * 120.0 / impulse;
			fy = fy * 120.0 / impulse;
		}

		// We apply the impact vector a second time
		this->vx -= fx / m1;
		this->vy -= fy / m1;
		u->vx += fx / m2;
		u->vy += fy / m2;
	}

	void output(const Move& move) {
		double a = angle + move.getAngle();

	    if (a >= 360.0) {
	        a = a - 360.0;
	    } else if (a < 0.0) {
	        a += 360.0;
	    }

	    // Look for a point corresponding to the angle we want
	    // Multiply to limit rounding errors
	    a = a * M_PI / 180.0;
	    double px = this->x + cos(a) * 100000.0;
	    double py = this->y + sin(a) * 100000.0;

	    if (move.getThrust() == -1) {
	    	this->activeShield() ;
	        std::cout << round(px) << " " << round(py) << " SHIELD" << std::endl;
	    }
	    else if (move.getThrust() == 650) {
	        this->useBoost();
	    	std::cout << round(px) << " " << round(py) << " BOOST" << std::endl;
	    }
		else {
			std::cout << round(px) << " " << round(py) << " " << move.getThrust() << std::endl;
	    }
	}

	Point moveToPoint(const Move& move) {
		double a = angle + move.getAngle();
		 if (a >= 360.0) {
			 a = a - 360.0;
		 } else if (a < 0.0) {
			 a += 360.0;
		 }
		 // Look for a point corresponding to the angle we want
		 // Multiply to limit rounding errors
		 a = a * M_PI / 180.0;
		 double px = round(this->x + cos(a) * 100000.0);
		 double py = round(this->y + sin(a) * 100000.0);

		return Point(px,py) ;
	}

	int getChecked() const {
		return checked;
	}

	double score(const Point* dest, const int max_checked) {
		if(this->getChecked() == max_checked) {
			return (max_checked +1) * 100 - this->getFinishTime() ;
		}
		return this->getChecked() * 100 - this->distanceSq(*dest) / MAX_DIST - abs(this->diffSpeedAngle(*dest)) / 180.0  ;

		/*return this->getChecked() * 100 - this->distanceSq(*dest) / MAX_DIST
				- std::min(abs(this->diffSpeedAngle(*dest)), abs(this->diffSpeedAngle(*nDest))) / 180.0 + speedSq / 100000000.0; */
	}

	bool getBoost() const {
		return hasBoost;
	}

	void useBoost() {
		this->hasBoost = false;
	}

protected:
	double angle ;
	int nextCheckpointId ;
	int timeout ;
	int shield ;
	int checked ;
	bool hasBoost ;
	double finishTime ;
};

class Checkpoint: public Unit {
public:
	Checkpoint(double x = -1., double y = -1.):
		Unit(x, y, 0, 0, 200.0)
	{//200 = 600 - 400
	}
};

class Game {
private:
	static int laps ;
	static int checkpointCount ;
	static int totalChecked ;

public:
	int turn ;
	static Checkpoint checkpoints[8] ;
	Pod pods[4] ;
	Game& operator=(Game const& g) ;

	Game() : turn(0) { }

	Game(int laps,int checkpointCount): turn(0) {
        Game::laps = laps ;
        Game::checkpointCount = checkpointCount ;
        Game::totalChecked = laps * checkpointCount ;

		for (int i = 0; i < checkpointCount; ++i) {
	    	int checkpointX;
	        int checkpointY;
	        std::cin >> checkpointX >> checkpointY; std::cin.ignore();
	        checkpoints[i] = Checkpoint(checkpointX, checkpointY) ;
	    }
	}

	Game(const Game& g) {
		turn = g.turn ;
		std::copy(std::begin(g.pods), std::end(g.pods), std::begin(this->pods)) ;
	}

	void readGame() {
		turn ++;

        for (int i = 0; i < 4; i++) {
            int x; // x position of your pod
            int y; // y position of your pod
            int vx; // x speed of your pod
            int vy; // y speed of your pod
            int angle; // angle of your pod
            int nextCheckPointId; // next check point id of your pod
            std::cin >> x >> y >> vx >> vy >> angle >> nextCheckPointId; std::cin.ignore();

            //warning : if you play p2, input angle is not -1 for opponent angle
            if(turn == 1) {
    			//if pod is not defined ; only on the first turn
                pods[i].updatePod(x,y,vx,vy,angle,nextCheckPointId) ;

            	//on first turn, pod angle is calculated based on first direction
    			//to avoid having a different behavior on first turn, we set the angle to next CP
    			double startAngle = pods[i].computeAngle(checkpoints[pods[i].getNextCheckpointId()]) ;
    			pods[i].setAngle(startAngle) ;
            }
            else {
            	if(nextCheckPointId != pods[i].getNextCheckpointId()) pods[i].checkedCP(this->nNextCP(pods + i)) ;
            	else pods[i].decTimeout() ;
            	pods[i].updatePod(x, y, vx, vy, angle, nextCheckPointId);
            	pods[i].decShield() ;
            }
        }
	}

	static const int getCheckpointCount() {
		return checkpointCount;
	}

	static const int getLaps() {
		return laps;
	}

	static int nNextCP(Pod* pod) {
		return pod->getNextCheckpointId() == checkpointCount - 1 ? 0 : pod->getNextCheckpointId() + 1 ;
	}

private:
	void simulateMovement() {
	    // This tracks the time during the turn. The goal is to reach 1.0
	    double t = 0.0;

	    while (t < 1.0) {
	        Collision* firstCollision = nullptr ;

	        // We look for all the collisions that are going to occur during the turn
	        for (int i = 0; i < 4 ; ++i) {
	            // Collision with another pod?
	            for (int j = i + 1 ; j < 4 ; ++j) {
	            	Collision* col = pods[i].collision(&pods[j]) ;

	                //we ignore the collision at t = 0
	                //even with them, I can't get an exact prediction
	                if(col != nullptr && col->t == 0) continue ;

	                // If the collision occurs earlier than the one we currently have we keep it
	                if (col != nullptr && col->t + t < 1.0 && (firstCollision == nullptr || col->t < firstCollision->t)) {
	                    firstCollision = col ;
	                }
	                else {
	                	delete col ;
	                }
	            }
	        }

            // Collision with checkpoints
	        //we use another loop for collision with CP, since it has no impact on Unit movements
	        for (int i = 0; i < 4 ; i++) {
	        	Collision* col = pods[i].collision(checkpoints + pods[i].getNextCheckpointId());

				// If the collision happens earlier than the current one play it
				if (col != nullptr && col->t + t < 1.0 && (firstCollision == nullptr || col->t < firstCollision->t)) {
					pods[i].checkedCP(this->nNextCP(pods+i)) ;
					//on sauvegarde le temps de franchissement du dernier CP
					if(pods[i].getChecked() == Game::totalChecked)
							pods[i].addFinishTime(col->t + t);
				}
				delete col ;
	        }

	        if (firstCollision == nullptr) {
	            // No collision, we can move the pods until the end of the turn
	            for (int i = 0; i < 4 ; ++i) {
	                pods[i].move(1.0 - t);
	            }

	            // End of the turn
	            t = 1.0;
	        } else {
	            // Move the pods to reach the time `t` of the collision
	            for (int i = 0; i < 4 ; ++i) {
	                pods[i].move(firstCollision->t);
	            }

	            // Play out the collision
            	static_cast<Pod*>(firstCollision->b)->bounce(static_cast<Pod*>(firstCollision->a)) ;

            	//this is needed to handle collision at t=0
            	//it adds a small value to avoid overlapping objects
            	//if(firstCollision->t == 0.0) Game::correctPositions(firstCollision) ;

	            t += firstCollision->t;
	        }
	    delete firstCollision ;
	    }

	    for (int i = 0; i < 4 ; ++i) {
	        pods[i].endTurn();
	        //on incrémente le temps mis pour franchir la ligne d'arrivée de 1
	        //la valeur de finish time vaut donc LENGTH toujours, sauf pour le dernier CP
	        if(pods[i].getChecked() < Game::totalChecked) pods[i].addFinishTime(1.0) ;
	    }

	    ++turn ;
	}

	static void correctPositions(Collision* col) {
		//specific test for collision at t=0
		//this happens when the distance between units is < to r+r
		//this is physically impossible but happens because of rounding errors
		double diff = col->a->distance(*col->b) - (col->a->getR() + col->b->getR()) ;
		//if(diff <= 0) {
			std::cerr << "Error before : " << diff << std::endl ;

			if(col->a->x > col->b->x) {
				col->a->x = col->a->x + EPSILON ;
				col->b->x = col->b->x - EPSILON ;
			}
			else {
				col->a->x = col->a->x - EPSILON ;
				col->b->x = col->b->x + EPSILON ;
			}

			if(col->a->y > col->b->y) {
				col->a->y = col->a->y + EPSILON ;
				col->b->y = col->b->y - EPSILON ;
			}
			else {
				col->a->y = col->a->y - EPSILON ;
				col->b->y = col->b->y + EPSILON ;
			}
			std::cerr << "Error after : " << col->a->distance(*col->b) - (col->a->getR() + col->b->getR()) << std::endl ;
		//}
	}

public:
	void simulateNextTurn(const Move& move0, const Move& move1, const Move& move2, const Move& move3) {
		//on arrête la simulation si au moins l'un des pods est arrivé
		/*const int max_checked = Game::laps * Game::checkpointCount ;
		if(pods[0].getChecked() == max_checked || pods[1].getChecked() == max_checked ||
				pods[2].getChecked() == max_checked || pods[3].getChecked() == max_checked) return ;*/

		this->pods[0].rotate(move0.getAngle());
		this->pods[0].boost(move0.getThrust()) ;

		this->pods[1].rotate(move1.getAngle());
		this->pods[1].boost(move1.getThrust()) ;

		this->pods[2].rotate(move2.getAngle());
		this->pods[2].boost(move2.getThrust()) ;

		this->pods[3].rotate(move3.getAngle());
		this->pods[3].boost(move3.getThrust()) ;

		this->simulateMovement();
	}

	void simulateNextTurn_ex(const Move& move0, const Move& move1, const Move& move2, const Move& move3) {
		this->pods[0].rotate(pods[0].moveToPoint(move0));
		this->pods[0].boost(move0.getThrust()) ;

		this->pods[1].rotate(pods[1].moveToPoint(move1));
		this->pods[1].boost(move1.getThrust()) ;

		this->pods[2].rotate(pods[2].moveToPoint(move2));
		this->pods[2].boost(move2.getThrust()) ;

		this->pods[3].rotate(pods[3].moveToPoint(move3));
		this->pods[3].boost(move3.getThrust()) ;

		this->simulateMovement();

        std::cerr << "x : " << pods[0].x << " y : " << pods[0].y << std::endl ;
        std::cerr << "vx : " << pods[0].getVx()  * 0.85 << " vy : " << pods[0].getVy()  * 0.85 << std::endl ;
       	std::cerr << "x : " << pods[1].x << " y : " << pods[1].y << std::endl ;
       	std::cerr << "vx : " << pods[1].getVx()  * 0.85  << " vy : " << pods[1].getVy() * 0.85  << std::endl ;
	}

	static double scorePod(Pod* pod){
		return pod->score(checkpoints + pod->getNextCheckpointId(), Game::totalChecked) ;
	}

	double evalGame(bool isOpp) {
		int fpi = 0, spi = 1, ofpi = 2, ospi =3 ;

		if(isOpp){
			fpi = 2;
			spi = 3 ;
			ofpi = 0 ;
			ospi = 1 ;
		}

		double score0 = scorePod(pods + fpi) ;
		double score1 = scorePod(pods + spi) ;
		double score2 = scorePod(pods + ofpi) ;
		double score3 = scorePod(pods + ospi) ;

		return score0 + score1 - std::max(score2, score3) ;
	}
};

Game& Game::operator=(Game const& g){
	turn = g.turn ;
	std::copy(std::begin(g.pods), std::end(g.pods), std::begin(this->pods)) ;
	return *this ;
}

//static variables initialization
int Game::laps = -1 ;
int Game::checkpointCount = -1 ;
int Game::totalChecked = -1 ;
Checkpoint Game::checkpoints[8] ;

class IA {
public:
	Game* game ;
	std::chrono::time_point<std::chrono::steady_clock> begin ;

	//virtual Solution* computeSolution() = 0 ;
	virtual ~IA() { } ;

	IA(Game* game, int max_time = DEFAULT_MAX_TIME, bool isOpp = false) :
		game(game), MAX_TIME(max_time), isOpp(isOpp) {
		begin = std::chrono::steady_clock::now();
	}

	void resetTimer() {
		begin = std::chrono::steady_clock::now();
	}

	double getElapsedTime() {
		return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - begin).count();
	}

	bool timeOut() {
		return getElapsedTime() >= MAX_TIME ;
	}

	void scoreSolution(Solution* s, Solution* os) {
		Game simulation = Game(*game) ;

		for(int i = 0 ; i < LENGTH ; ++i) {
			simulation.simulateNextTurn(s->solution[i], s->solution[i+LENGTH], os->solution[i], os->solution[i+LENGTH]) ;
		}
		if(!isOpp) s->setScore(simulation.evalGame(false)) ;
		else os->setScore(simulation.evalGame(true)) ;
	}

protected :
	const int MAX_TIME ;
	const bool isOpp ;
};

class SimpleIA : public IA {
public:
	SimpleIA(Game* game, bool isOpp=false) : IA(game,0,isOpp) {
		if(!isOpp) {
			myPods = game->pods ;
			oppPods = game->pods + 2 ;
		}
		else {
			myPods = game->pods +2 ;
			oppPods = game -> pods ;
		}
	}

	~SimpleIA() {

	}
protected:
	//advance to next checkpoint
	Move computeAMove(Pod* pod) {
		Unit* target = game->checkpoints + pod->getNextCheckpointId() ;

		double speed = sqrt(pod->getVx() * pod->getVx() + pod->getVy() * pod->getVy()) ;
		double thrust = 200.0 ;
		double n = pod->distance(*target) * 0.85 / speed ;
		double angle = pod->diffAngle(*target) ;
		double correction = 0. ;

		//new target : next checkpoint
		Unit* newTarget = game->checkpoints + Game::nNextCP(pod) ;
		double newAngle = pod->diffAngle(*newTarget) ;

		//if we are close to the checkpoint, we try to go to the next one
		//this allows to turn before reaching checkpoint
		if(n < 3.0 && abs(newAngle) > 45 && pod->getChecked() < game->getCheckpointCount() * game->getLaps() - 1) {
			Pod clone = Pod(*pod) ;
			for(int j = 0 ; j < n ; ++j) {
				clone.rotate(*newTarget) ;
				//boost 0 is useless
				//clone.boost(0.) ;
				Collision* col = clone.collision(target) ;
				if(col != nullptr) {
					target = newTarget ;
					angle = newAngle ;
					thrust = 0 ;
					delete col ;
					goto end ;
				}
				else {
					clone.move(1.0) ;
					clone.endTurn() ;
				}
			}
		}

		//we try to correct the angle
		//correction is calculated to counteract the angle between direction and speed
		if(abs(angle) < 18.0 && speed != 0 && thrust != 0) {
			correction = 2 * pod->diffSpeedAngle(*target) * thrust / speed ;
		}
		else if(abs(angle) > 90.0) {
			thrust = 0 ;
		}

		if(abs(angle) > 36 && abs(pod->diffSpeedAngle(*target)) > 45) thrust = 0. ;

		end:
		angle = angle + correction ;
		if(angle > 18.0) angle = 18.0 ;
		else if(angle < -18.0) angle = -18.0 ;
		return Move(angle, thrust) ;
	}

	//block opponent
	Move computeBMove(Pod* pod, Pod* opp) {
		//pod is our blocker, opp is the opponent to block
		double dist2opp = pod->distanceSq(*opp) ;
		Point oppDest = Point(opp->getX() + opp->getVx(), opp->getY() + opp->getVy()) ;

		Point* target = nullptr ;
		double thrust = 0 ;

		//on fonce sur l'adv à proximité
		if(dist2opp < 5000 * 5000) {
			target = &oppDest ;
			thrust = 200. ;
		}
		//si a proximité du CP de l'adv, on attend
		else if(pod->distanceSq(game->checkpoints[opp->getNextCheckpointId()]) < 2000 * 2000) {
			target = &oppDest ;
			thrust = 0. ;
		}
		else if(pod->distanceSq(game->checkpoints[Game::nNextCP(opp)]) < 2000 * 2000) {
			target = &oppDest ;
			thrust = 0. ;
		}
		else if(pod->distanceSq(game->checkpoints[opp->getNextCheckpointId()]) < opp->distanceSq(game->checkpoints[opp->getNextCheckpointId()])){
			target = game->checkpoints + opp->getNextCheckpointId() ;
			thrust = 100 ;
		}
		else {
			target = game->checkpoints + Game::nNextCP(opp) ;
			thrust = 100 ;
		}

		return Move(pod->diffAngle(*target),thrust) ;
	}

	bool checkCollision(Pod* p, const Move& m, Pod* o) {
		Pod pod = Pod(*p) ;
		pod.rotate(m.getAngle()) ;
		pod.boost(m.getThrust()) ;

		Pod opp = Pod(*o) ;
		opp.rotate(game->checkpoints[opp.getNextCheckpointId()]) ;
		opp.boost(200.);

		Collision* col = pod.collision(&opp) ;
		if(col != nullptr) {
			delete col ;
			return true ;
		}
		return false ;
	}

public:
	Move* computeMoves() {
		Move* moves = new Move[2] ;

		//best opponent
		Pod* opp = std::max(oppPods,oppPods+1,
				[this](Pod* a,Pod* b) {
					return Game::scorePod(a) < Game::scorePod(b) ;
				}
		);

		//bool firstTurn = myPods[0].getChecked() < game->getCheckpointCount() &&
		//		myPods[1].getChecked() < game->getCheckpointCount() ;

		double score0 = Game::scorePod(myPods) ;
		double score1 = Game::scorePod(myPods + 1) ;

		/*if(firstTurn) {
			moves[0] = computeAMove(myPods) ;
			moves[1] = computeAMove(myPods + 1) ;
		}
		else if(score0 > score1) {*/
		if(score0 > score1) {
			moves[0] = computeAMove(myPods) ;
			moves[1] = computeBMove(myPods + 1, opp) ;
		}
		else {
			moves[0] = computeBMove(myPods, opp) ;
			moves[1] = computeAMove(myPods + 1) ;
		}

		if(checkCollision(myPods, moves[0], oppPods) || checkCollision(myPods, moves[0], oppPods +1)) moves[0].setThrust(-1.) ;
		if(checkCollision(myPods+1, moves[1], oppPods) || checkCollision(myPods+1, moves[1], oppPods+1)) moves[1].setThrust(-1.) ;

		return moves ;
	}

	template <class IA>
	static Solution computeSolution(Game* g, bool isOpp){
		Solution s ;
		Game simulation(*g) ;
		IA ia(&simulation, isOpp) ;

		//simulation des tours consécutifs pour construction de la solution
		if(!isOpp) {
			for(int i = 0 ; i < LENGTH ; ++i) {
				Move* m(ia.computeMoves()) ;
				s.solution[i] = m[0] ;
				s.solution[i+LENGTH] = m[1] ;
				simulation.simulateNextTurn(m[0], m[1], Move(0,0), Move(0,0)) ;
				delete[] m ;
			}
		}
		else {
			for(int i = 0 ; i < LENGTH ; ++i) {
				Move* m(ia.computeMoves()) ;
				s.solution[i] = m[0] ;
				s.solution[i+LENGTH] = m[1] ;
				simulation.simulateNextTurn(Move(0,0), Move(0,0), m[0], m[1]) ;
				delete[] m ;
			}
		}
		return s ;
	}

protected:
	Pod *myPods ;
	Pod *oppPods ;
};

//IA : the two pods try to race as fast as possible
class RaceIA : public SimpleIA {
public:
	RaceIA(Game* game, bool isOpp=false) : SimpleIA(game, isOpp) {

	}
	~RaceIA() { }

	Move* computeMoves() {
		Move* mtab = new Move[2] ;
		mtab[0] = computeAMove(myPods) ;
		mtab[1] = computeAMove(myPods + 1) ;

		if(checkCollision(myPods, mtab[0], oppPods) || checkCollision(myPods, mtab[0], oppPods +1)) mtab[0].setThrust(-1.) ;
		if(checkCollision(myPods+1, mtab[1], oppPods) || checkCollision(myPods+1, mtab[1], oppPods+1)) mtab[1].setThrust(-1.) ;

		return mtab ;
	}
};

class BlockIA : public SimpleIA {
public:
	BlockIA(Game* game, bool isOpp=false) : SimpleIA(game, isOpp) {

	}
	~BlockIA() { }

	Move* computeMoves() {
		//best opponent
		Pod* opp = std::max(oppPods,oppPods+1,
				[this](Pod* a,Pod* b) {
					return Game::scorePod(a) < Game::scorePod(b) ;
				}
		);

		Move* mtab = new Move[2] ;
		mtab[0] = computeBMove(myPods, opp) ;
		mtab[1] = computeBMove(myPods + 1, opp) ;

		if(checkCollision(myPods, mtab[0], oppPods) || checkCollision(myPods, mtab[0], oppPods +1)) mtab[0].setThrust(-1.) ;
		if(checkCollision(myPods+1, mtab[1], oppPods) || checkCollision(myPods+1, mtab[1], oppPods+1)) mtab[1].setThrust(-1.) ;

		return mtab ;
	}
};

class SAIA : public IA {
private:
	const int NUM_ITERATION = 100 ; //number of iteration at each T°
	const double alpha = 0.97 ; //temperature reduction at each iteration
	const int INITIAL_TEMP = 10 ; //initial T°
	const bool KEEP_BEST = false ;

	long total_iterations ;
	bool hasBestSolution ;

	Solution bestSolution ;

	static bool acceptance(double oldValue, double newValue, double temp) {
		Random& r = Random::getInstance();
		double deltaE = newValue - oldValue ;

		if(deltaE >= 0) return true ;
		else {
			//boltzmann acceptance
			if(r.nextFloat() < exp(deltaE/temp)) return true ;
			else return false ;
		}
	}

public:
	SAIA(Game* game, int max_time = DEFAULT_MAX_TIME, bool isOpp=false) : IA(game, max_time, isOpp),
			total_iterations(0), hasBestSolution(false), bestSolution() {
	}

	~SAIA() { }

	long getTotalIterations() const {
		return total_iterations ;
	}

	Solution computeSolution(Solution& opponentSolution) {
		double temp = INITIAL_TEMP ;

		//initial solution building & scoring
		if(!KEEP_BEST || !hasBestSolution) {
			bestSolution = SimpleIA::computeSolution<RaceIA>(game, isOpp) ;
			if(!isOpp) this->scoreSolution(&bestSolution, &opponentSolution);
			else this->scoreSolution(&opponentSolution,&bestSolution);
			std::cerr << "Initial solution at temp " << temp << " with score " << bestSolution.getScore() << std::endl;
		}
		else {
			bestSolution.shiftLeft() ;
			if(!isOpp) this->scoreSolution(&bestSolution, &opponentSolution);
			else this->scoreSolution(&opponentSolution,&bestSolution);
			std::cerr << "Reusing best solution at temp " << temp << " with score " << bestSolution.getScore() << std::endl;
		}

		hasBestSolution = true ;
		Solution currentSolution(bestSolution) ;

		while(!timeOut()) {
			for(int i = 0 ; i < NUM_ITERATION && !timeOut() ; ++i) {

				//creation of new solution
				currentSolution.mutate() ;
				++total_iterations ;

				//score the new solution
				if(!isOpp) this->scoreSolution(&currentSolution, &opponentSolution);
				else this->scoreSolution(&opponentSolution, &currentSolution);

				//is new solution accepted ?
				if(acceptance(currentSolution.getSavedScore(), currentSolution.getScore(), temp)) {
					if(currentSolution.getScore() > bestSolution.getScore()) {
						bestSolution = currentSolution ;
					}
				}
				else {
					currentSolution.restore() ;
				}
			}
			temp = temp * alpha ;
		}
		std::cerr << "Best solution at temp " << temp << " with score " << bestSolution.getScore() << std::endl;
		return bestSolution ;
	}

	Move* computeMoves(Solution& opponentSolution) {
		Solution s = this->computeSolution(opponentSolution) ;
		Move* moves = new Move[2] ;
		moves[0] = s.solution[0] ;
		moves[1] = s.solution[LENGTH] ;
		return moves ;
	}

};

int main()
{
	std::ios::sync_with_stdio(false);
    int laps;
    std::cin >> laps; std::cin.ignore();
    int checkpointCount;
    std::cin >> checkpointCount; std::cin.ignore();

    Game game(laps,checkpointCount) ;
    SAIA opponentIA(&game,45,true) ;
    SAIA ia(&game,100) ;
    Move* moves = nullptr ;

    // game loop
    while (1) {
    	game.readGame();
    	opponentIA.resetTimer() ;
    	Solution basic = SimpleIA::computeSolution<SimpleIA>(&game, false) ;
    	Solution opponentSolution = opponentIA.computeSolution(basic) ;

    	ia.resetTimer();
    	//Solution opponentSolution = SimpleIA::computeSolution<SimpleIA>(&game, true) ;
    	moves = ia.computeMoves(opponentSolution) ;

    	std::cerr << "Time : " << ia.getElapsedTime() << std::endl ;
    	std::cerr << "Average # of iterations : " << ia.getTotalIterations() / game.turn << std::endl ;

//     	Game clone(game) ;
//      clone.simulateNextTurn_ex(moves[0],moves[1],Move(0,80),Move(0,80));

    	game.pods[0].output(moves[0]);
    	game.pods[1].output(moves[1]);
    	delete[] moves ;
    }

    delete[] moves ;
}
