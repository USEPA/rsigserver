
#ifndef ASSERTIONS_H
#define ASSERTIONS_H

/******************************************************************************
PURPOSE: Assertions.h - Defines assertion macros for optionally verifying
         routine preconditions, postconditions and class invariants -
         as inspired by Eiffel's 'Design By Contract'.
         http://www.elj.com/eiffel/dbc/

NOTES:   This implementation is based on (uses) the ISO/ANSI Standard C/C++
         Pre-processor (cpp).

         Macros defined include:

           assert_static(c)             Assert condition at compile time.
                                        Note this does not generate code at
                                        runtime and therefore is unaffected by
                                        the various assertion level options
                                        defined below.

           PRE0(c), PRE02(c1,c2), ...   Assert precondition(s)
           POST0(c), POST02(c1,c2), ... Assert postcondition(s)
           PRE(c),  PRE2(c1,c2),  ...   Assert invariant and precondition(s)
           POST(c), POST2(c1,c2), ...   Assert invariant and postcondition(s)

         There are three levels of assertion checking and they are
         specified on the command line or in the Makefile:

           (default)            - Enables  all assertion checking.
           -DNO_ASSERTIONS      - Disables all assertion checking.
           -DPRECONDITIONS_ONLY - Enables PRE(c) checking only.

         Specifying -DNO_INVARIANT during pre-processing causes the PRE() and
         POST() macros to NOT also assert the class invariant, otherwise, by
         default, they invoke a function assumed to be in local (module)
         scope with the following signature:

         If C++:

           virtual bool invariant() const;

         If C:

           int invariant( self );

         Use -DINVARIANT_SELF=THIS to rename the first argument of each
         module call to 'THIS' instead of the default, 'self'.

         In C++ constructors or code in headers or any non-ADT code,
         use PRE0, PRE02, POST0, ..., macros to avoid invariant evaluation.

         Additional assertion checking:

           CHECK(c), CHECK2(c1,c2),       Assert arbitrary condition(s)
                                          Useful to insert in the middle of
                                          a long routine to verify
                                          intermediate results.

           CHECKING(s)                    Includes arbitrary statements, s,
                                          only when CHECK macros are
                                          enabled. This is useful for
                                          including/removing variables and
                                          routines whose only use occurs
                                          within CHECK/POST macros.

         Similaraly, there are two levels of debugging:

           (default)            - Disables DEBUG(s) statements.
           -DDEBUGGING          - Enables  DEBUG(s) statements.

           DEBUG(s)                       Add optional debug statements
                                          independent of assertions.

         Assertions support macros include:

           OLD(variable)                  For referring to previous values.
           REMEMBER(type,variable)        Required if OLD() is used...
           REMEMBER_F(type,function_name) Same as above but for functions.


         Macros for boolean expressions. The functional forms e.g.,

           IMPLIES_ELSE( found,
                         AND2( index > 0, items( index ) == sought_item ),
                         index == -1 )

         have been shown to be more readable than the infix-operator form:

           found && ( index > 0 && items( index ) == sought_item )
           || ! found && index == -1


           IMPLIES(p,c)                   Logical implication.
           IMPLIES_ELSE(p,c1,c2)          p -> c1 or !p ->c2
           AND2(c1,c2), AND3(c1,c2,c3)... Logical and.
           OR2(c1,c2),  OR3(c1,c2,c3)...  Logical or.
           XOR2(c1,c2), XOR3(c1,c2,c3)... Logical exclusive-or.
                                          XORn(c1,c2,...) is 1 iff exactly
                                          one of the conditions is true.
                                          Not short-circuited.
           IN3(x,a,b),  IN4(x,a,b,c)...   Set membership: x in {a,b}?
           IN_RANGE(x,low,high)           x in [low,high]?
           IS_BOOL2(a,b), ...             Each in {0, 1}?
           NON_ZERO2(a,b), ...            Each != 0?
           IS_ZERO2(a,b), ...             Each == 0?
           GT_ZERO2(a,b), ...             Each >  0?
           GE_ZERO2(a,b), ...             Each >= 0?


         If any asserted condition evaluates to false at runtime then
         the program will display a message such as:

           Assertion failed: PRE: index >= 0, file Test.c, line 44

         And then abort() producing a core file used for debugging.
         For example, to debug a bad call to a routine (violated
         pre-condition):

           dbx Test
           dbx> up 4
           (this is the assertion line)
           dbx> up
           (this is the bad call line)

         Checking preconditions or postconditions imply checking the
         invariant also. If there is no precondition or postcondition then
         use PRE( true ) and POST( true ) so that the invariant is still
         checked. See search() below for an example of this.

         The invariant should be a non-pure virtual function
         (one with a body that at least returns true) so that it can be
         called by decendents. See invariant() below.

         First C++ Example:

         $ ls
         Assertions.h        assertions_test.C
         $ cat assertions_test.C

         // assertions_test.C - demonstrate Assertions.h macros.

         #include <iostream> // For std::cout, std::endl.

         #include <Assertions.h> // For PRE*(), POST*(), CHECK*(), IN_RANGE().

         class Counter {
         private:
           int value;   // Current value of counter.
           int maximum; // Maximum value allowed.
         public:

           ~Counter() {
             PRE( true ); // Use PRE( true ) so that invariant() is evaluated.
             value = maximum = 0;
           }

           Counter( int maximum_value = 10 )
             : value( 1 ), maximum( maximum_value ) {
             POST2( count() == 1, limit() == maximum_value ); // 'anded' args.
           }

           int limit() const {
             PRE( true );
             const int result = maximum;
             POST( result > 0 );
             return result;
           }

           int count() const {
             PRE( true );
             const int result = value;
             POST( IN_RANGE( result, 1, limit() ) ); // IN_RANGE is inclusive.
             return result;
           }

           void increment() {
             PRE( count() < limit() );
             REMEMBER_F( int, count ); // declares 'int count_old_ = count();'
             ++value;
             POST( count() == OLD( count ) + 1 );
           }

           void decrement() {
             PRE( count() > 1 );
             REMEMBER_F( int, count );
             --value;
             POST( count() == OLD( count ) - 1 );
           }

           // Class invariant: AND2(limit() > 0, IN_RANGE(count(), 1, limit()))
           bool invariant() const {
             const bool result = AND2( maximum > 0,
                                       IN_RANGE( value, 1, maximum ) );
             return result;
           }
         };

         int main() {
           assert_static( sizeof (long long) == 8 );
           assert_static( sizeof (double) == 8 );
           Counter c;
           std::cout << "c.count() = " << c.count() << std::endl;
           std::cout << "c.limit() = " << c.limit() << std::endl;
           c.increment();
           CHECK( c.count() == 2 );
           c.increment();
           CHECK( c.count() == 3 );
           c.decrement();
           CHECK( c.count() == 2 );
           std::cout << "c.count() = " << c.count() << std::endl;
           c.decrement();
           CHECK( c.count() == 1 );
           c.decrement(); // Defect: decrement when count() is 1.
           std::cout << "c.count() = " << c.count() << std::endl;
           return 0;
         }

         $ CC -LANG:std -I. -g -o assertions_test assertions_test.C
         $ assertions_test
         c.count() = 1
         c.limit() = 10
         c.count() = 2
         Assertion failed: PRE: count() > 1, file assertions_test.C, line 51
         Abort (core dumped)
         $ dbx assertions_test
         dbx version 7.3.1 68542_Oct26 MR Oct 26 2000 17:50:34
         Core from signal SIGABRT: Abort (see abort(3c))
         (dbx) up 4
         Counter::decrement:  51  PRE( count() > 1 );
         (dbx) up
         ::main:  77  c.decrement(); // Defect: decrement when count() is 1.
         (dbx) quit
         $ \rm core
         $ CC -LANG:std -I. -DNO_ASSERTIONS -O \
              -o assertions_test assertions_test.C
         $ assertions_test
         c.count() = 1
         c.limit() = 10
         c.count() = 2
         c.count() = 0


         Second C++ example with inheritance:

           // Array10.h:

           class Array10 : public Array, public FixedSize {
           private:
             int n;
             double a[ 10 ];
           public:
             virtual bool invariant() const; // This is called by PRE/POST().
             bool insert( int i, double x ); // Insert x at a[ i ].
             int search(  double x ); // Return index of x in a[] or -1.
             // ...
           };

           // Array10.C:

           #include <Assertions.h>
           #include <Array10.h>

           bool Array10::invariant() const {
             // Should be written to call each of the immediate parent's
             // invariant() also (no variance):

             return AND3( Array::invariant(), FixedSize::invariant(),
                          IN_RANGE( n, 0, 10 ) );
                          // The IN_RANGE() part is this class's invariant.
           }

           void Array10::insert( int i, double x ) {
             PRE2( 0 <= i, i <= n );
             REMEMBER( int, n ); // REMEMBER() here if using OLD() in POST().

             // ...
             a[ i ] = x;
             // ...

             POST2( a[ i ] == x, n == OLD( n ) + 1 ); // Also invariant.
           }

           int Array10::search( double x ) {
             PRE( true ); // Must check the invariant.
             int result = -1;
             size_t i;

             for ( i = 0; i < n; ++i ) if ( a[ i ] == x ) break;

             if ( i < n ) result = i;

             POST( IMPLIES( result != -1, a[ result ] == x ) );
             return result;
           }

           The CHECK() macro can be used to assert a condition without
           implying the class invariant - e.g., inside a helper routine.
           The DEBUG() macro is used to optionally insert any statements.

           Examples:

           int X::helper( int i ) {
             CHECK( i > 0 );
             DEBUG( cout << "X::helper( int i = " << i << " )" << endl; )
             // ...
             CHECK( result > 0 );
             DEBUG( cout << "X::helper() returning " << result << endl; )
             return result;
           }

           When compiling with C, the PRE*() macros declare and initialize
           a variable 'unused_assert_' that is otherwise not used.
           This is done to enable the PRE*() macros to be placed at the
           beginning of the routine (just like C++). For example:

           int half_fast_routine( int value ) {
             PRE( value >= 0 );
             const int result = value >> 1;
             POST2( result >= 0, result + result <= value );
             return result;
           }

           With or without assertions enabled, the semi-colon following PRE()
           and before the declaration of result requires that PRE() be some
           kind of declaration. When assertions are enabled it declares a
           local unused variable unused_assert_.
           When assertions are disabled, it declares 'extern int errno'
           (a global variable in the Standard C Library) which does not cause
           any code-generation. Groovy.

HISTORY: 1994-02-01 plessel@computer.org
STATUS:  reviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#ifndef NO_ASSERTIONS
#include <assert.h> /* For function void __assert(). */
#endif

/*================================= MACROS ==================================*/

/*
 * Define a compile-time assertion that only accepts compile-time constant
 * expressions and if false results in a compile-time error and if true
 * generates no instructions or warnings.
 *
 * Example usage: assert_static( sizeof (long long) == 8 );
 */

#define assert_static(a) extern int unused_assert_static_[ (a) ? 1 : -1 ]


/*
 * Define a custom assert-like macro that includes a prefix tag that identifies
 * the assertion type (PRE, POST, CHECK).
 * This macro calls the "Standard" C Library function void __assert().
 */

#ifndef NO_ASSERTIONS

#if defined( __MINGW32__ ) || defined( __MINGW64__ )
#define ASSERT2_(prefix,expression) \
((void) ((expression) ? 0 \
:(_assert(#prefix": "#expression,__FILE__,__LINE__), 0)))
#elif defined( _WIN32 )
#define ASSERT2_(prefix,expression) \
((void) ((expression) ? 0 \
:(_assert(__FILE__,__LINE__,#prefix": "#expression), 0)))
#elif defined(__CYGWIN__)
#define ASSERT2_(prefix,expression) \
((void) ((expression) ? 0 \
:(__assert(__FILE__,__LINE__,#prefix": "#expression), 0)))
#else
#define ASSERT2_(prefix,expression) \
((void) ((expression) ? 0 \
:(__assert(#prefix": "#expression,__FILE__,__LINE__), 0)))
#endif

#define PRE_ASSERT_(expression) ASSERT2_(PRE,expression)
#define POST_ASSERT_(expression) ASSERT2_(POST,expression)
#define CHECK_ASSERT_(expression) ASSERT2_(CHECK,expression)
#else
#define ASSERT2_(unused1,unused2)
#define PRE_ASSERT_(unused)
#define POST_ASSERT_(unused)
#define CHECK_ASSERT_(unused)
#endif


#ifdef NO_ASSERTIONS

/******************************* NO_ASSERTIONS *******************************/

/* Disable all checking. */

#define CHECK(c)
#define CHECK2(c1,c2)
#define CHECK3(c1,c2,c3)
#define CHECK4(c1,c2,c3,c4)
#define CHECK5(c1,c2,c3,c4,c5)
#define CHECK6(c1,c2,c3,c4,c5,c6)
#define CHECK7(c1,c2,c3,c4,c5,c6,c7)
#define CHECK8(c1,c2,c3,c4,c5,c6,c7,c8)
#define CHECK9(c1,c2,c3,c4,c5,c6,c7,c8,c9)
#define CHECK10(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)
#define CHECK11(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11)
#define CHECK12(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12)
#define CHECK13(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13)
#define CHECK14(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14)
#define CHECK15(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15)
#define CHECK16(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16)
#define CHECK17(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17)
#define CHECK18(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18)

#define CHECK19(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19)

#define CHECK20(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20)

#define CHECK21(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21)

#define CHECK22(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22)

#define CHECK23(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23)

#define CHECK24(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24)

#define CHECK25(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25)

#define CHECK26(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26)

#define CHECK27(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27)

#define CHECK28(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28)

#define CHECK29(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29)

#define CHECK30(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30)

/* An external variable declaration is needed in C so pre can be first line. */

#ifdef __cplusplus
#define NOP_DECLARATION_
#else
#define NOP_DECLARATION_ extern int errno
#endif

#define PRE(c1) NOP_DECLARATION_
#define PRE2(c1,c2) NOP_DECLARATION_
#define PRE3(c1,c2,c3) NOP_DECLARATION_
#define PRE4(c1,c2,c3,c4) NOP_DECLARATION_
#define PRE5(c1,c2,c3,c4,c5) NOP_DECLARATION_
#define PRE6(c1,c2,c3,c4,c5,c6) NOP_DECLARATION_
#define PRE7(c1,c2,c3,c4,c5,c6,c7) NOP_DECLARATION_
#define PRE8(c1,c2,c3,c4,c5,c6,c7,c8) NOP_DECLARATION_
#define PRE9(c1,c2,c3,c4,c5,c6,c7,c8,c9) NOP_DECLARATION_
#define PRE10(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10) NOP_DECLARATION_
#define PRE11(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11) NOP_DECLARATION_
#define PRE12(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12) NOP_DECLARATION_
#define PRE13(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13) NOP_DECLARATION_
#define PRE14(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14) NOP_DECLARATION_

#define PRE15(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15) \
 NOP_DECLARATION_

#define PRE16(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16) \
 NOP_DECLARATION_

#define PRE17(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17) \
 NOP_DECLARATION_

#define PRE18(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18) \
 NOP_DECLARATION_

#define PRE19(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19) NOP_DECLARATION_

#define PRE20(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20) NOP_DECLARATION_

#define PRE21(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21) NOP_DECLARATION_

#define PRE22(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22) NOP_DECLARATION_

#define PRE23(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23) NOP_DECLARATION_

#define PRE24(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24) NOP_DECLARATION_

#define PRE25(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25) NOP_DECLARATION_

#define PRE26(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26) NOP_DECLARATION_

#define PRE27(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27) NOP_DECLARATION_

#define PRE28(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28) NOP_DECLARATION_

#define PRE29(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29) NOP_DECLARATION_

#define PRE30(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30) NOP_DECLARATION_

#define PRE0 PRE
#define PRE02 PRE2
#define PRE03 PRE3
#define PRE04 PRE4
#define PRE05 PRE5
#define PRE06 PRE6
#define PRE07 PRE7
#define PRE08 PRE8
#define PRE09 PRE9
#define PRE010 PRE10
#define PRE011 PRE11
#define PRE012 PRE12
#define PRE013 PRE13
#define PRE014 PRE14
#define PRE015 PRE15
#define PRE016 PRE16
#define PRE017 PRE17
#define PRE018 PRE18
#define PRE019 PRE19
#define PRE020 PRE20
#define PRE021 PRE21
#define PRE022 PRE22
#define PRE023 PRE23
#define PRE024 PRE24
#define PRE025 PRE25
#define PRE026 PRE26
#define PRE027 PRE27
#define PRE028 PRE28
#define PRE029 PRE29
#define PRE030 PRE30

#define POST CHECK
#define POST2 CHECK2
#define POST3 CHECK3
#define POST4 CHECK4
#define POST5 CHECK5
#define POST6 CHECK6
#define POST7 CHECK7
#define POST8 CHECK8
#define POST9 CHECK9
#define POST10 CHECK10
#define POST11 CHECK11
#define POST12 CHECK12
#define POST13 CHECK13
#define POST14 CHECK14
#define POST15 CHECK15
#define POST16 CHECK16
#define POST17 CHECK17
#define POST18 CHECK18
#define POST19 CHECK19
#define POST20 CHECK20
#define POST21 CHECK21
#define POST22 CHECK22
#define POST23 CHECK23
#define POST24 CHECK24
#define POST25 CHECK25
#define POST26 CHECK26
#define POST27 CHECK27
#define POST28 CHECK28
#define POST29 CHECK29
#define POST30 CHECK30

#define POST0 CHECK
#define POST02 CHECK2
#define POST03 CHECK3
#define POST04 CHECK4
#define POST05 CHECK5
#define POST06 CHECK6
#define POST07 CHECK7
#define POST08 CHECK8
#define POST09 CHECK9
#define POST010 CHECK10
#define POST011 CHECK11
#define POST012 CHECK12
#define POST013 CHECK13
#define POST014 CHECK14
#define POST015 CHECK15
#define POST016 CHECK16
#define POST017 CHECK17
#define POST018 CHECK18
#define POST019 CHECK19
#define POST020 CHECK20
#define POST021 CHECK21
#define POST022 CHECK22
#define POST023 CHECK23
#define POST024 CHECK24
#define POST025 CHECK25
#define POST026 CHECK26
#define POST027 CHECK27
#define POST028 CHECK28
#define POST029 CHECK29
#define POST030 CHECK30

#define CHECKING(s)

#define OLD(variable)
#define REMEMBER(type,variable) NOP_DECLARATION_
#define REMEMBER_F(type,function_name) NOP_DECLARATION_


#elif defined(PRECONDITIONS_ONLY)

/***************************** PRECONDITIONS_ONLY ****************************/

/* Enable pre-condition checking only. */

#ifdef __cplusplus

#define PRE(c) PRE_ASSERT_(c)

#define PRE2(c1,c2) PRE_ASSERT_(c1),PRE_ASSERT_(c2)

#define PRE3(c1,c2,c3) PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3)

#define PRE4(c1,c2,c3,c4) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4)

#define PRE5(c1,c2,c3,c4,c5) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),PRE_ASSERT_(c5)

#define PRE6(c1,c2,c3,c4,c5,c6) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6)

#define PRE7(c1,c2,c3,c4,c5,c6,c7) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7)

#define PRE8(c1,c2,c3,c4,c5,c6,c7,c8) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8)

#define PRE9(c1,c2,c3,c4,c5,c6,c7,c8,c9) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9)

#define PRE10(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),PRE_ASSERT_(c10)

#define PRE11(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11)

#define PRE12(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),
PRE_ASSERT_(c12)

#define PRE13(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),\
PRE_ASSERT_(c12),PRE_ASSERT_(c13)

#define PRE14(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),\
PRE_ASSERT_(c12),PRE_ASSERT_(c13),PRE_ASSERT_(c14)

#define PRE15(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15)

#define PRE16(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16)

#define PRE17(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17)

#define PRE18(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18)

#define PRE19(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19)

#define PRE20(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20)

#define PRE21(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21)

#define PRE22(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22)

#define PRE23(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23)

#define PRE24(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24)

#define PRE25(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25)

#define PRE26(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25),PRE_ASSERT_(c26)

#define PRE27(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25),PRE_ASSERT_(c26),PRE_ASSERT_(c27)

#define PRE28(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25),PRE_ASSERT_(c26),PRE_ASSERT_(c27),PRE_ASSERT_(c28)

#define PRE29(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25),PRE_ASSERT_(c26),PRE_ASSERT_(c27),PRE_ASSERT_(c28),\
PRE_ASSERT_(c29)

#define PRE30(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25),PRE_ASSERT_(c26),PRE_ASSERT_(c27),PRE_ASSERT_(c28),\
PRE_ASSERT_(c29),PRE_ASSERT_(c30)

#else /* ! __cplusplus */

/* C versions declare a local variable so they can be the first statement. */

#define PRE(c) const int unused_assert_ = (PRE_ASSERT_(c),&unused_assert_!=0)

#define PRE2(c1,c2) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),&unused_assert_!=0)

#define PRE3(c1,c2,c3) \
const int unused_assert_ = \
(PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),&unused_assert_!=0)

#define PRE4(c1,c2,c3,c4) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),&unused_assert_!=0)

#define PRE5(c1,c2,c3,c4,c5) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),&unused_assert_!=0)

#define PRE6(c1,c2,c3,c4,c5,c6) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),\
PRE_ASSERT_(c6),&unused_assert_!=0)

#define PRE7(c1,c2,c3,c4,c5,c6,c7) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
&unused_assert_!=0)

#define PRE8(c1,c2,c3,c4,c5,c6,c7,c8) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),&unused_assert_!=0)

#define PRE9(c1,c2,c3,c4,c5,c6,c7,c8,c9) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),&unused_assert_!=0)

#define PRE10(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),PRE_ASSERT_(c10),&unused_assert_!=0)

#define PRE11(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),\
&unused_assert_!=0)

#define PRE12(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),\
PRE_ASSERT_(c12),&unused_assert_!=0)

#define PRE13(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),\
PRE_ASSERT_(c12),PRE_ASSERT_(c13),&unused_assert_!=0)

#define PRE14(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),\
PRE_ASSERT_(c12),PRE_ASSERT_(c13),PRE_ASSERT_(c14),&unused_assert_!=0)

#define PRE15(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),\
PRE_ASSERT_(c12),PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),\
&unused_assert_!=0)

#define PRE16(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),\
PRE_ASSERT_(c12),PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),\
PRE_ASSERT_(c16),&unused_assert_!=0)

#define PRE17(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),\
PRE_ASSERT_(c12),PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),\
PRE_ASSERT_(c16),PRE_ASSERT_(c17),&unused_assert_!=0)

#define PRE18(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),\
PRE_ASSERT_(c12),PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),\
PRE_ASSERT_(c16),PRE_ASSERT_(c17),PRE_ASSERT_(c18),&unused_assert_!=0)

#define PRE19(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),\
PRE_ASSERT_(c12),PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),\
PRE_ASSERT_(c16),PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),\
&unused_assert_!=0)

#define PRE20(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),\
PRE_ASSERT_(c12),PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),\
PRE_ASSERT_(c16),PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),\
PRE_ASSERT_(c20),&unused_assert_!=0)

#define PRE21(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),\
PRE_ASSERT_(c12),PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),\
PRE_ASSERT_(c16),PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),\
PRE_ASSERT_(c20),PRE_ASSERT_(c21),&unused_assert_!=0)

#define PRE22(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),\
PRE_ASSERT_(c12),PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),\
PRE_ASSERT_(c16),PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),\
PRE_ASSERT_(c20),PRE_ASSERT_(c21),PRE_ASSERT_(c22),&unused_assert_!=0)

#define PRE23(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),\
PRE_ASSERT_(c12),PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),\
PRE_ASSERT_(c16),PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),\
PRE_ASSERT_(c20),PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),\
&unused_assert_!=0)

#define PRE24(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),\
PRE_ASSERT_(c12),PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),\
PRE_ASSERT_(c16),PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),\
PRE_ASSERT_(c20),PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),\
PRE_ASSERT_(c24),&unused_assert_!=0)

#define PRE25(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),\
PRE_ASSERT_(c12),PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),\
PRE_ASSERT_(c16),PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),\
PRE_ASSERT_(c20),PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),\
PRE_ASSERT_(c24),PRE_ASSERT_(c25),&unused_assert_!=0)

#define PRE26(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),\
PRE_ASSERT_(c12),PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),\
PRE_ASSERT_(c16),PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),\
PRE_ASSERT_(c20),PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),\
PRE_ASSERT_(c24),PRE_ASSERT_(c25),PRE_ASSERT_(c26),&unused_assert_!=0)

#define PRE27(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),\
PRE_ASSERT_(c12),PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),\
PRE_ASSERT_(c16),PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),\
PRE_ASSERT_(c20),PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),\
PRE_ASSERT_(c24),PRE_ASSERT_(c25),PRE_ASSERT_(c26),PRE_ASSERT_(c27),\
&unused_assert_!=0)

#define PRE28(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),\
PRE_ASSERT_(c12),PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),\
PRE_ASSERT_(c16),PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),\
PRE_ASSERT_(c20),PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),\
PRE_ASSERT_(c24),PRE_ASSERT_(c25),PRE_ASSERT_(c26),PRE_ASSERT_(c27),\
PRE_ASSERT_(c28),&unused_assert_!=0)

#define PRE29(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),\
PRE_ASSERT_(c12),PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),\
PRE_ASSERT_(c16),PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),\
PRE_ASSERT_(c20),PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),\
PRE_ASSERT_(c24),PRE_ASSERT_(c25),PRE_ASSERT_(c26),PRE_ASSERT_(c27),\
PRE_ASSERT_(c28),PRE_ASSERT_(c29),&unused_assert_!=0)

#define PRE30(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),\
PRE_ASSERT_(c4),PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),\
PRE_ASSERT_(c8),PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),\
PRE_ASSERT_(c12),PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),\
PRE_ASSERT_(c16),PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),\
PRE_ASSERT_(c20),PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),\
PRE_ASSERT_(c24),PRE_ASSERT_(c25),PRE_ASSERT_(c26),PRE_ASSERT_(c27),\
PRE_ASSERT_(c28),PRE_ASSERT_(c29),PRE_ASSERT_(c30),&unused_assert_!=0)


#endif /*!  __cplusplus */


#define CHECK(c)
#define CHECK2(c1,c2)
#define CHECK3(c1,c2,c3)
#define CHECK4(c1,c2,c3,c4)
#define CHECK5(c1,c2,c3,c4,c5)
#define CHECK6(c1,c2,c3,c4,c5,c6)
#define CHECK7(c1,c2,c3,c4,c5,c6,c7)
#define CHECK8(c1,c2,c3,c4,c5,c6,c7,c8)
#define CHECK9(c1,c2,c3,c4,c5,c6,c7,c8,c9)
#define CHECK10(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)
#define CHECK11(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11)
#define CHECK12(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12)
#define CHECK13(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13)
#define CHECK14(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14)
#define CHECK15(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15)
#define CHECK16(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16)
#define CHECK17(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17)
#define CHECK18(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18)

#define CHECK19(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19)

#define CHECK20(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20)

#define CHECK21(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21)

#define CHECK22(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22)

#define CHECK23(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23)

#define CHECK24(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24)

#define CHECK25(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25)

#define CHECK26(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26)

#define CHECK27(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27)

#define CHECK28(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28)

#define CHECK29(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29)

#define CHECK30(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30)

#define PRE0 PRE
#define PRE02 PRE2
#define PRE03 PRE3
#define PRE04 PRE4
#define PRE05 PRE5
#define PRE06 PRE6
#define PRE07 PRE7
#define PRE08 PRE8
#define PRE09 PRE9
#define PRE010 PRE10
#define PRE011 PRE11
#define PRE012 PRE12
#define PRE013 PRE13
#define PRE014 PRE14
#define PRE015 PRE15
#define PRE016 PRE16
#define PRE017 PRE17
#define PRE018 PRE18
#define PRE019 PRE19
#define PRE020 PRE20
#define PRE021 PRE21
#define PRE022 PRE22
#define PRE023 PRE23
#define PRE024 PRE24
#define PRE025 PRE25
#define PRE026 PRE26
#define PRE027 PRE27
#define PRE028 PRE28
#define PRE029 PRE29
#define PRE030 PRE30

#define POST CHECK
#define POST2 CHECK2
#define POST3 CHECK3
#define POST4 CHECK4
#define POST5 CHECK5
#define POST6 CHECK6
#define POST7 CHECK7
#define POST8 CHECK8
#define POST9 CHECK9
#define POST10 CHECK10
#define POST11 CHECK11
#define POST12 CHECK12
#define POST13 CHECK13
#define POST14 CHECK14
#define POST15 CHECK15
#define POST16 CHECK16
#define POST17 CHECK17
#define POST18 CHECK18
#define POST19 CHECK19
#define POST20 CHECK20
#define POST21 CHECK21
#define POST22 CHECK22
#define POST23 CHECK23
#define POST24 CHECK24
#define POST25 CHECK25
#define POST26 CHECK26
#define POST27 CHECK27
#define POST28 CHECK28
#define POST29 CHECK29
#define POST30 CHECK30

#define POST0 CHECK
#define POST02 CHECK2
#define POST03 CHECK3
#define POST04 CHECK4
#define POST05 CHECK5
#define POST06 CHECK6
#define POST07 CHECK7
#define POST08 CHECK8
#define POST09 CHECK9
#define POST010 CHECK10
#define POST011 CHECK11
#define POST012 CHECK12
#define POST013 CHECK13
#define POST014 CHECK14
#define POST015 CHECK15
#define POST016 CHECK16
#define POST017 CHECK17
#define POST018 CHECK18
#define POST019 CHECK19
#define POST020 CHECK20
#define POST021 CHECK21
#define POST022 CHECK22
#define POST023 CHECK23
#define POST024 CHECK24
#define POST025 CHECK25
#define POST026 CHECK26
#define POST027 CHECK27
#define POST028 CHECK28
#define POST029 CHECK29
#define POST030 CHECK30

#define CHECKING(s)

#define NOP_DECLARATION_ extern int errno

#define OLD(variable)
#define REMEMBER(type,variable) NOP_DECLARATION_
#define REMEMBER_F(type,function_name) NOP_DECLARATION_


#else /* ! PRECONDITIONS_ONLY */

/******************************* ALL ASSERTIONS ******************************/

/* Enable all checking. */

#define CHECK(c) CHECK_ASSERT_(c)

#define CHECK2(c1,c2) CHECK_ASSERT_(c1),CHECK_ASSERT_(c2)

#define CHECK3(c1,c2,c3) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3)

#define CHECK4(c1,c2,c3,c4) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4)

#define CHECK5(c1,c2,c3,c4,c5) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5)

#define CHECK6(c1,c2,c3,c4,c5,c6) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6)

#define CHECK7(c1,c2,c3,c4,c5,c6,c7) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7)

#define CHECK8(c1,c2,c3,c4,c5,c6,c7,c8) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8)

#define CHECK9(c1,c2,c3,c4,c5,c6,c7,c8,c9) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9)

#define CHECK10(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9),CHECK_ASSERT_(c10)

#define CHECK11(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9),CHECK_ASSERT_(c10),CHECK_ASSERT_(c11)

#define CHECK12(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9),CHECK_ASSERT_(c10),CHECK_ASSERT_(c11),CHECK_ASSERT_(c12)

#define CHECK13(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9),CHECK_ASSERT_(c10),CHECK_ASSERT_(c11),CHECK_ASSERT_(c12),\
CHECK_ASSERT_(c13)

#define CHECK14(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9),CHECK_ASSERT_(c10),CHECK_ASSERT_(c11),CHECK_ASSERT_(c12),\
CHECK_ASSERT_(c13),CHECK_ASSERT_(c14)

#define CHECK15(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9),CHECK_ASSERT_(c10),CHECK_ASSERT_(c11),CHECK_ASSERT_(c12),\
CHECK_ASSERT_(c13),CHECK_ASSERT_(c14),CHECK_ASSERT_(c15)

#define CHECK16(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9),CHECK_ASSERT_(c10),CHECK_ASSERT_(c11),CHECK_ASSERT_(c12),\
CHECK_ASSERT_(c13),CHECK_ASSERT_(c14),CHECK_ASSERT_(c15),CHECK_ASSERT_(c16)

#define CHECK17(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9),CHECK_ASSERT_(c10),CHECK_ASSERT_(c11),CHECK_ASSERT_(c12),\
CHECK_ASSERT_(c13),CHECK_ASSERT_(c14),CHECK_ASSERT_(c15),CHECK_ASSERT_(c16),\
CHECK_ASSERT_(c17)

#define CHECK18(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9),CHECK_ASSERT_(c10),CHECK_ASSERT_(c11),CHECK_ASSERT_(c12),\
CHECK_ASSERT_(c13),CHECK_ASSERT_(c14),CHECK_ASSERT_(c15),CHECK_ASSERT_(c16),\
CHECK_ASSERT_(c17),CHECK_ASSERT_(c18)

#define CHECK19(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9),CHECK_ASSERT_(c10),CHECK_ASSERT_(c11),CHECK_ASSERT_(c12),\
CHECK_ASSERT_(c13),CHECK_ASSERT_(c14),CHECK_ASSERT_(c15),CHECK_ASSERT_(c16),\
CHECK_ASSERT_(c17),CHECK_ASSERT_(c18),CHECK_ASSERT_(c19)

#define CHECK20(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9),CHECK_ASSERT_(c10),CHECK_ASSERT_(c11),CHECK_ASSERT_(c12),\
CHECK_ASSERT_(c13),CHECK_ASSERT_(c14),CHECK_ASSERT_(c15),CHECK_ASSERT_(c16),\
CHECK_ASSERT_(c17),CHECK_ASSERT_(c18),CHECK_ASSERT_(c19),CHECK_ASSERT_(c20)

#define CHECK21(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9),CHECK_ASSERT_(c10),CHECK_ASSERT_(c11),CHECK_ASSERT_(c12),\
CHECK_ASSERT_(c13),CHECK_ASSERT_(c14),CHECK_ASSERT_(c15),CHECK_ASSERT_(c16),\
CHECK_ASSERT_(c17),CHECK_ASSERT_(c18),CHECK_ASSERT_(c19),CHECK_ASSERT_(c20),\
CHECK_ASSERT_(c21)

#define CHECK22(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9),CHECK_ASSERT_(c10),CHECK_ASSERT_(c11),CHECK_ASSERT_(c12),\
CHECK_ASSERT_(c13),CHECK_ASSERT_(c14),CHECK_ASSERT_(c15),CHECK_ASSERT_(c16),\
CHECK_ASSERT_(c17),CHECK_ASSERT_(c18),CHECK_ASSERT_(c19),CHECK_ASSERT_(c20),\
CHECK_ASSERT_(c21),CHECK_ASSERT_(c22)

#define CHECK23(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9),CHECK_ASSERT_(c10),CHECK_ASSERT_(c11),CHECK_ASSERT_(c12),\
CHECK_ASSERT_(c13),CHECK_ASSERT_(c14),CHECK_ASSERT_(c15),CHECK_ASSERT_(c16),\
CHECK_ASSERT_(c17),CHECK_ASSERT_(c18),CHECK_ASSERT_(c19),CHECK_ASSERT_(c20),\
CHECK_ASSERT_(c21),CHECK_ASSERT_(c22),CHECK_ASSERT_(c23)

#define CHECK24(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9),CHECK_ASSERT_(c10),CHECK_ASSERT_(c11),CHECK_ASSERT_(c12),\
CHECK_ASSERT_(c13),CHECK_ASSERT_(c14),CHECK_ASSERT_(c15),CHECK_ASSERT_(c16),\
CHECK_ASSERT_(c17),CHECK_ASSERT_(c18),CHECK_ASSERT_(c19),CHECK_ASSERT_(c20),\
CHECK_ASSERT_(c21),CHECK_ASSERT_(c22),CHECK_ASSERT_(c23),CHECK_ASSERT_(c24)

#define CHECK25(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9),CHECK_ASSERT_(c10),CHECK_ASSERT_(c11),CHECK_ASSERT_(c12),\
CHECK_ASSERT_(c13),CHECK_ASSERT_(c14),CHECK_ASSERT_(c15),CHECK_ASSERT_(c16),\
CHECK_ASSERT_(c17),CHECK_ASSERT_(c18),CHECK_ASSERT_(c19),CHECK_ASSERT_(c20),\
CHECK_ASSERT_(c21),CHECK_ASSERT_(c22),CHECK_ASSERT_(c23),CHECK_ASSERT_(c24),\
CHECK_ASSERT_(c25)

#define CHECK26(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9),CHECK_ASSERT_(c10),CHECK_ASSERT_(c11),CHECK_ASSERT_(c12),\
CHECK_ASSERT_(c13),CHECK_ASSERT_(c14),CHECK_ASSERT_(c15),CHECK_ASSERT_(c16),\
CHECK_ASSERT_(c17),CHECK_ASSERT_(c18),CHECK_ASSERT_(c19),CHECK_ASSERT_(c20),\
CHECK_ASSERT_(c21),CHECK_ASSERT_(c22),CHECK_ASSERT_(c23),CHECK_ASSERT_(c24),\
CHECK_ASSERT_(c25),CHECK_ASSERT_(c26)

#define CHECK27(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9),CHECK_ASSERT_(c10),CHECK_ASSERT_(c11),CHECK_ASSERT_(c12),\
CHECK_ASSERT_(c13),CHECK_ASSERT_(c14),CHECK_ASSERT_(c15),CHECK_ASSERT_(c16),\
CHECK_ASSERT_(c17),CHECK_ASSERT_(c18),CHECK_ASSERT_(c19),CHECK_ASSERT_(c20),\
CHECK_ASSERT_(c21),CHECK_ASSERT_(c22),CHECK_ASSERT_(c23),CHECK_ASSERT_(c24),\
CHECK_ASSERT_(c25),CHECK_ASSERT_(c26),CHECK_ASSERT_(c27)

#define CHECK28(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9),CHECK_ASSERT_(c10),CHECK_ASSERT_(c11),CHECK_ASSERT_(c12),\
CHECK_ASSERT_(c13),CHECK_ASSERT_(c14),CHECK_ASSERT_(c15),CHECK_ASSERT_(c16),\
CHECK_ASSERT_(c17),CHECK_ASSERT_(c18),CHECK_ASSERT_(c19),CHECK_ASSERT_(c20),\
CHECK_ASSERT_(c21),CHECK_ASSERT_(c22),CHECK_ASSERT_(c23),CHECK_ASSERT_(c24),\
CHECK_ASSERT_(c25),CHECK_ASSERT_(c26),CHECK_ASSERT_(c27),CHECK_ASSERT_(c28)

#define CHECK29(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9),CHECK_ASSERT_(c10),CHECK_ASSERT_(c11),CHECK_ASSERT_(c12),\
CHECK_ASSERT_(c13),CHECK_ASSERT_(c14),CHECK_ASSERT_(c15),CHECK_ASSERT_(c16),\
CHECK_ASSERT_(c17),CHECK_ASSERT_(c18),CHECK_ASSERT_(c19),CHECK_ASSERT_(c20),\
CHECK_ASSERT_(c21),CHECK_ASSERT_(c22),CHECK_ASSERT_(c23),CHECK_ASSERT_(c24),\
CHECK_ASSERT_(c25),CHECK_ASSERT_(c26),CHECK_ASSERT_(c27),CHECK_ASSERT_(c28),\
CHECK_ASSERT_(c29)

#define CHECK30(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9),CHECK_ASSERT_(c10),CHECK_ASSERT_(c11),CHECK_ASSERT_(c12),\
CHECK_ASSERT_(c13),CHECK_ASSERT_(c14),CHECK_ASSERT_(c15),CHECK_ASSERT_(c16),\
CHECK_ASSERT_(c17),CHECK_ASSERT_(c18),CHECK_ASSERT_(c19),CHECK_ASSERT_(c20),\
CHECK_ASSERT_(c21),CHECK_ASSERT_(c22),CHECK_ASSERT_(c23),CHECK_ASSERT_(c24),\
CHECK_ASSERT_(c25),CHECK_ASSERT_(c26),CHECK_ASSERT_(c27),CHECK_ASSERT_(c28),\
CHECK_ASSERT_(c29),CHECK_ASSERT_(c30)



#define POST0(c) POST_ASSERT_(c)

#define POST02(c1,c2) POST_ASSERT_(c1),POST_ASSERT_(c2)

#define POST03(c1,c2,c3) POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3)

#define POST04(c1,c2,c3,c4) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4)

#define POST05(c1,c2,c3,c4,c5) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5)

#define POST06(c1,c2,c3,c4,c5,c6) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6)

#define POST07(c1,c2,c3,c4,c5,c6,c7) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7)

#define POST08(c1,c2,c3,c4,c5,c6,c7,c8) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8)

#define POST09(c1,c2,c3,c4,c5,c6,c7,c8,c9) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9)

#define POST010(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9),POST_ASSERT_(c10)

#define POST011(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9),POST_ASSERT_(c10),POST_ASSERT_(c11)

#define POST012(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9),POST_ASSERT_(c10),POST_ASSERT_(c11),POST_ASSERT_(c12)

#define POST013(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9),POST_ASSERT_(c10),POST_ASSERT_(c11),POST_ASSERT_(c12),\
POST_ASSERT_(c13)

#define POST014(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9),POST_ASSERT_(c10),POST_ASSERT_(c11),POST_ASSERT_(c12),\
POST_ASSERT_(c13),POST_ASSERT_(c14)

#define POST015(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9),POST_ASSERT_(c10),POST_ASSERT_(c11),POST_ASSERT_(c12),\
POST_ASSERT_(c13),POST_ASSERT_(c14),POST_ASSERT_(c15)

#define POST016(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9),POST_ASSERT_(c10),POST_ASSERT_(c11),POST_ASSERT_(c12),\
POST_ASSERT_(c13),POST_ASSERT_(c14),POST_ASSERT_(c15),POST_ASSERT_(c16)

#define POST017(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9),POST_ASSERT_(c10),POST_ASSERT_(c11),POST_ASSERT_(c12),\
POST_ASSERT_(c13),POST_ASSERT_(c14),POST_ASSERT_(c15),POST_ASSERT_(c16),\
POST_ASSERT_(c17)

#define POST018(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9),POST_ASSERT_(c10),POST_ASSERT_(c11),POST_ASSERT_(c12),\
POST_ASSERT_(c13),POST_ASSERT_(c14),POST_ASSERT_(c15),POST_ASSERT_(c16),\
POST_ASSERT_(c17),POST_ASSERT_(c18)

#define POST019(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9),POST_ASSERT_(c10),POST_ASSERT_(c11),POST_ASSERT_(c12),\
POST_ASSERT_(c13),POST_ASSERT_(c14),POST_ASSERT_(c15),POST_ASSERT_(c16),\
POST_ASSERT_(c17),POST_ASSERT_(c18),POST_ASSERT_(c19)

#define POST020(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9),POST_ASSERT_(c10),POST_ASSERT_(c11),POST_ASSERT_(c12),\
POST_ASSERT_(c13),POST_ASSERT_(c14),POST_ASSERT_(c15),POST_ASSERT_(c16),\
POST_ASSERT_(c17),POST_ASSERT_(c18),POST_ASSERT_(c19),POST_ASSERT_(c20)

#define POST021(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9),POST_ASSERT_(c10),POST_ASSERT_(c11),POST_ASSERT_(c12),\
POST_ASSERT_(c13),POST_ASSERT_(c14),POST_ASSERT_(c15),POST_ASSERT_(c16),\
POST_ASSERT_(c17),POST_ASSERT_(c18),POST_ASSERT_(c19),POST_ASSERT_(c20),\
POST_ASSERT_(c21)

#define POST022(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9),POST_ASSERT_(c10),POST_ASSERT_(c11),POST_ASSERT_(c12),\
POST_ASSERT_(c13),POST_ASSERT_(c14),POST_ASSERT_(c15),POST_ASSERT_(c16),\
POST_ASSERT_(c17),POST_ASSERT_(c18),POST_ASSERT_(c19),POST_ASSERT_(c20),\
POST_ASSERT_(c21),POST_ASSERT_(c22)

#define POST023(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9),POST_ASSERT_(c10),POST_ASSERT_(c11),POST_ASSERT_(c12),\
POST_ASSERT_(c13),POST_ASSERT_(c14),POST_ASSERT_(c15),POST_ASSERT_(c16),\
POST_ASSERT_(c17),POST_ASSERT_(c18),POST_ASSERT_(c19),POST_ASSERT_(c20),\
POST_ASSERT_(c21),POST_ASSERT_(c22),POST_ASSERT_(c23)

#define POST024(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9),POST_ASSERT_(c10),POST_ASSERT_(c11),POST_ASSERT_(c12),\
POST_ASSERT_(c13),POST_ASSERT_(c14),POST_ASSERT_(c15),POST_ASSERT_(c16),\
POST_ASSERT_(c17),POST_ASSERT_(c18),POST_ASSERT_(c19),POST_ASSERT_(c20),\
POST_ASSERT_(c21),POST_ASSERT_(c22),POST_ASSERT_(c23),POST_ASSERT_(c24)

#define POST025(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9),POST_ASSERT_(c10),POST_ASSERT_(c11),POST_ASSERT_(c12),\
POST_ASSERT_(c13),POST_ASSERT_(c14),POST_ASSERT_(c15),POST_ASSERT_(c16),\
POST_ASSERT_(c17),POST_ASSERT_(c18),POST_ASSERT_(c19),POST_ASSERT_(c20),\
POST_ASSERT_(c21),POST_ASSERT_(c22),POST_ASSERT_(c23),POST_ASSERT_(c24),\
POST_ASSERT_(c25)

#define POST026(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9),POST_ASSERT_(c10),POST_ASSERT_(c11),POST_ASSERT_(c12),\
POST_ASSERT_(c13),POST_ASSERT_(c14),POST_ASSERT_(c15),POST_ASSERT_(c16),\
POST_ASSERT_(c17),POST_ASSERT_(c18),POST_ASSERT_(c19),POST_ASSERT_(c20),\
POST_ASSERT_(c21),POST_ASSERT_(c22),POST_ASSERT_(c23),POST_ASSERT_(c24),\
POST_ASSERT_(c25),POST_ASSERT_(c26)

#define POST027(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9),POST_ASSERT_(c10),POST_ASSERT_(c11),POST_ASSERT_(c12),\
POST_ASSERT_(c13),POST_ASSERT_(c14),POST_ASSERT_(c15),POST_ASSERT_(c16),\
POST_ASSERT_(c17),POST_ASSERT_(c18),POST_ASSERT_(c19),POST_ASSERT_(c20),\
POST_ASSERT_(c21),POST_ASSERT_(c22),POST_ASSERT_(c23),POST_ASSERT_(c24),\
POST_ASSERT_(c25),POST_ASSERT_(c26),POST_ASSERT_(c27)

#define POST028(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9),POST_ASSERT_(c10),POST_ASSERT_(c11),POST_ASSERT_(c12),\
POST_ASSERT_(c13),POST_ASSERT_(c14),POST_ASSERT_(c15),POST_ASSERT_(c16),\
POST_ASSERT_(c17),POST_ASSERT_(c18),POST_ASSERT_(c19),POST_ASSERT_(c20),\
POST_ASSERT_(c21),POST_ASSERT_(c22),POST_ASSERT_(c23),POST_ASSERT_(c24),\
POST_ASSERT_(c25),POST_ASSERT_(c26),POST_ASSERT_(c27),POST_ASSERT_(c28)

#define POST029(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9),POST_ASSERT_(c10),POST_ASSERT_(c11),POST_ASSERT_(c12),\
POST_ASSERT_(c13),POST_ASSERT_(c14),POST_ASSERT_(c15),POST_ASSERT_(c16),\
POST_ASSERT_(c17),POST_ASSERT_(c18),POST_ASSERT_(c19),POST_ASSERT_(c20),\
POST_ASSERT_(c21),POST_ASSERT_(c22),POST_ASSERT_(c23),POST_ASSERT_(c24),\
POST_ASSERT_(c25),POST_ASSERT_(c26),POST_ASSERT_(c27),POST_ASSERT_(c28),\
POST_ASSERT_(c29)

#define POST030(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9),POST_ASSERT_(c10),POST_ASSERT_(c11),POST_ASSERT_(c12),\
POST_ASSERT_(c13),POST_ASSERT_(c14),POST_ASSERT_(c15),POST_ASSERT_(c16),\
POST_ASSERT_(c17),POST_ASSERT_(c18),POST_ASSERT_(c19),POST_ASSERT_(c20),\
POST_ASSERT_(c21),POST_ASSERT_(c22),POST_ASSERT_(c23),POST_ASSERT_(c24),\
POST_ASSERT_(c25),POST_ASSERT_(c26),POST_ASSERT_(c27),POST_ASSERT_(c28),\
POST_ASSERT_(c29),POST_ASSERT_(c30)



#ifdef __cplusplus

#ifndef NO_INVARIANT
#define EVALUATE_INVARIANT CHECK_ASSERT_(invariant())
#else
#define EVALUATE_INVARIANT ((void)(0))
#endif


#define PRE0(c) PRE_ASSERT_(c)

#define PRE02(c1,c2) PRE_ASSERT_(c1),PRE_ASSERT_(c2)

#define PRE03(c1,c2,c3) PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3)

#define PRE04(c1,c2,c3,c4) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4)

#define PRE05(c1,c2,c3,c4,c5) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5)

#define PRE06(c1,c2,c3,c4,c5,c6) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6)

#define PRE07(c1,c2,c3,c4,c5,c6,c7) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7)

#define PRE08(c1,c2,c3,c4,c5,c6,c7,c8) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8)

#define PRE09(c1,c2,c3,c4,c5,c6,c7,c8,c9) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9)

#define PRE010(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10)

#define PRE011(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11)

#define PRE012(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12)

#define PRE013(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13)

#define PRE014(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14)

#define PRE015(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15)

#define PRE016(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16)

#define PRE017(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17)

#define PRE018(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18)

#define PRE019(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19)

#define PRE020(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20)

#define PRE021(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21)

#define PRE022(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22)

#define PRE023(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23)

#define PRE024(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24)

#define PRE025(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25)

#define PRE026(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25),PRE_ASSERT_(c26)

#define PRE027(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25),PRE_ASSERT_(c26),PRE_ASSERT_(c27)

#define PRE028(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25),PRE_ASSERT_(c26),PRE_ASSERT_(c27),PRE_ASSERT_(c28)

#define PRE029(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25),PRE_ASSERT_(c26),PRE_ASSERT_(c27),PRE_ASSERT_(c28),\
PRE_ASSERT_(c29)

#define PRE030(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30) \
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25),PRE_ASSERT_(c26),PRE_ASSERT_(c27),PRE_ASSERT_(c28),\
PRE_ASSERT_(c29),PRE_ASSERT_(c30)


#define PRE(c) EVALUATE_INVARIANT,PRE0(c)

#define PRE2(c1,c2) EVALUATE_INVARIANT,PRE02(c1,c2)

#define PRE3(c1,c2,c3) EVALUATE_INVARIANT,PRE03(c1,c2,c3)

#define PRE4(c1,c2,c3,c4) EVALUATE_INVARIANT,PRE04(c1,c2,c3,c4)

#define PRE5(c1,c2,c3,c4,c5) EVALUATE_INVARIANT,PRE05(c1,c2,c3,c4,c5)

#define PRE6(c1,c2,c3,c4,c5,c6) EVALUATE_INVARIANT,PRE06(c1,c2,c3,c4,c5,c6)

#define PRE7(c1,c2,c3,c4,c5,c6,c7) \
EVALUATE_INVARIANT,PRE07(c1,c2,c3,c4,c5,c6,c7)

#define PRE8(c1,c2,c3,c4,c5,c6,c7,c8) \
EVALUATE_INVARIANT,PRE08(c1,c2,c3,c4,c5,c6,c7,c8)

#define PRE9(c1,c2,c3,c4,c5,c6,c7,c8,c9) \
EVALUATE_INVARIANT,PRE09(c1,c2,c3,c4,c5,c6,c7,c8,c9)

#define PRE10(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10) \
EVALUATE_INVARIANT,PRE010(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)

#define PRE11(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11) \
EVALUATE_INVARIANT,PRE011(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11)

#define PRE12(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12) \
EVALUATE_INVARIANT,PRE012(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12)

#define PRE13(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13) \
EVALUATE_INVARIANT,PRE013(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13)

#define PRE14(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14) \
EVALUATE_INVARIANT,PRE014(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14)

#define PRE15(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15) \
EVALUATE_INVARIANT,PRE015(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15)

#define PRE16(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16) \
EVALUATE_INVARIANT,PRE016(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,\
c16)

#define PRE17(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17) \
EVALUATE_INVARIANT,\
PRE017(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17)

#define PRE18(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,\
c18) \
EVALUATE_INVARIANT,\
PRE018(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18)

#define PRE19(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19) \
EVALUATE_INVARIANT,\
PRE019(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19)

#define PRE20(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20) \
EVALUATE_INVARIANT,\
PRE020(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20)

#define PRE21(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,\
c18,c19,c20,c21) \
EVALUATE_INVARIANT,\
PRE021(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,\
c21)

#define PRE22(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22) \
EVALUATE_INVARIANT,\
PRE022(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,\
c21,c22)

#define PRE23(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23) \
EVALUATE_INVARIANT,\
PRE023(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,\
c21,c22,c23)

#define PRE24(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24) \
EVALUATE_INVARIANT,\
PRE024(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,\
c21,c22,c23,c24)

#define PRE25(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25) \
EVALUATE_INVARIANT,\
PRE025(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,\
c21,c22,c23,c24,c25)

#define PRE26(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26) \
EVALUATE_INVARIANT,\
PRE026(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,\
c21,c22,c23,c24,c25,c26)

#define PRE27(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27) \
EVALUATE_INVARIANT,\
PRE027(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,\
c21,c22,c23,c24,c25,c26,c27)

#define PRE28(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28) \
EVALUATE_INVARIANT,\
PRE028(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,\
c21,c22,c23,c24,c25,c26,c27,c28)

#define PRE29(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29) \
EVALUATE_INVARIANT,\
PRE029(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,\
c21,c22,c23,c24,c25,c26,c27,c28,c29)

#define PRE30(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30) \
EVALUATE_INVARIANT,\
PRE030(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,\
c21,c22,c23,c24,c25,c26,c27,c28,c29,c30)

#else /* ! __cplusplus */


#ifndef NO_INVARIANT
#ifndef INVARIANT_SELF
#define INVARIANT_SELF self
#endif
#define EVALUATE_INVARIANT \
  CHECK_ASSERT_(INVARIANT_SELF?INVARIANT_SELF->invariant(INVARIANT_SELF):0)
#else
#define EVALUATE_INVARIANT ((void)(0))
#endif


#define PRE0(c) const int unused_assert_ = (PRE_ASSERT_(c),&unused_assert_!=0)

#define PRE02(c1,c2) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),&unused_assert_!=0)

#define PRE03(c1,c2,c3) \
const int unused_assert_ = \
(PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),&unused_assert_!=0)

#define PRE04(c1,c2,c3,c4) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
&unused_assert_!=0)

#define PRE05(c1,c2,c3,c4,c5) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),&unused_assert_!=0)

#define PRE06(c1,c2,c3,c4,c5,c6) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),&unused_assert_!=0)

#define PRE07(c1,c2,c3,c4,c5,c6,c7) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),&unused_assert_!=0)

#define PRE08(c1,c2,c3,c4,c5,c6,c7,c8) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
&unused_assert_!=0)

#define PRE09(c1,c2,c3,c4,c5,c6,c7,c8,c9) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),&unused_assert_!=0)

#define PRE010(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),&unused_assert_!=0)

#define PRE011(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),&unused_assert_!=0)

#define PRE012(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
&unused_assert_!=0)

#define PRE013(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),&unused_assert_!=0)

#define PRE014(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),&unused_assert_!=0)

#define PRE015(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),&unused_assert_!=0)

#define PRE016(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
&unused_assert_!=0)

#define PRE017(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),&unused_assert_!=0)

#define PRE018(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),&unused_assert_!=0)

#define PRE019(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),&unused_assert_!=0)

#define PRE020(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
&unused_assert_!=0)

#define PRE021(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),&unused_assert_!=0)

#define PRE022(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),&unused_assert_!=0)

#define PRE023(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),&unused_assert_!=0)

#define PRE024(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
&unused_assert_!=0)

#define PRE025(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25),&unused_assert_!=0)

#define PRE026(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25),PRE_ASSERT_(c26),&unused_assert_!=0)

#define PRE027(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25),PRE_ASSERT_(c26),PRE_ASSERT_(c27),&unused_assert_!=0)

#define PRE028(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25),PRE_ASSERT_(c26),PRE_ASSERT_(c27),PRE_ASSERT_(c28),\
&unused_assert_!=0)

#define PRE029(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25),PRE_ASSERT_(c26),PRE_ASSERT_(c27),PRE_ASSERT_(c28),\
PRE_ASSERT_(c29), &unused_assert_!=0)

#define PRE030(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25),PRE_ASSERT_(c26),PRE_ASSERT_(c27),PRE_ASSERT_(c28),\
PRE_ASSERT_(c29),PRE_ASSERT_(c30), &unused_assert_!=0)



#define PRE(c) const int unused_assert_ = \
(EVALUATE_INVARIANT,PRE_ASSERT_(c),&unused_assert_!=0)

#define PRE2(c1,c2) \
const int unused_assert_ = \
(EVALUATE_INVARIANT,PRE_ASSERT_(c1),PRE_ASSERT_(c2),&unused_assert_!=0)

#define PRE3(c1,c2,c3) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),&unused_assert_!=0)

#define PRE4(c1,c2,c3,c4) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
&unused_assert_!=0)

#define PRE5(c1,c2,c3,c4,c5) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),&unused_assert_!=0)

#define PRE6(c1,c2,c3,c4,c5,c6) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),&unused_assert_!=0)

#define PRE7(c1,c2,c3,c4,c5,c6,c7) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),&unused_assert_!=0)

#define PRE8(c1,c2,c3,c4,c5,c6,c7,c8) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
&unused_assert_!=0)

#define PRE9(c1,c2,c3,c4,c5,c6,c7,c8,c9) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),&unused_assert_!=0)

#define PRE10(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),&unused_assert_!=0)

#define PRE11(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),&unused_assert_!=0)

#define PRE12(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
&unused_assert_!=0)

#define PRE13(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),&unused_assert_!=0)

#define PRE14(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),&unused_assert_!=0)

#define PRE15(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),&unused_assert_!=0)

#define PRE16(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
&unused_assert_!=0)

#define PRE17(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),&unused_assert_!=0)

#define PRE18(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),&unused_assert_!=0)

#define PRE19(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),&unused_assert_!=0)

#define PRE20(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
&unused_assert_!=0)

#define PRE21(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),&unused_assert_!=0)

#define PRE22(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),&unused_assert_!=0)

#define PRE23(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),&unused_assert_!=0)

#define PRE24(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
&unused_assert_!=0)

#define PRE25(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25),&unused_assert_!=0)

#define PRE26(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25),PRE_ASSERT_(c26),&unused_assert_!=0)

#define PRE27(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25),PRE_ASSERT_(c26),PRE_ASSERT_(c27),&unused_assert_!=0)

#define PRE28(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25),PRE_ASSERT_(c26),PRE_ASSERT_(c27),PRE_ASSERT_(c28),\
&unused_assert_!=0)

#define PRE29(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25),PRE_ASSERT_(c26),PRE_ASSERT_(c27),PRE_ASSERT_(c28),\
PRE_ASSERT_(c29), &unused_assert_!=0)

#define PRE30(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30) \
const int unused_assert_ = (EVALUATE_INVARIANT,\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),PRE_ASSERT_(c10),PRE_ASSERT_(c11),PRE_ASSERT_(c12),\
PRE_ASSERT_(c13),PRE_ASSERT_(c14),PRE_ASSERT_(c15),PRE_ASSERT_(c16),\
PRE_ASSERT_(c17),PRE_ASSERT_(c18),PRE_ASSERT_(c19),PRE_ASSERT_(c20),\
PRE_ASSERT_(c21),PRE_ASSERT_(c22),PRE_ASSERT_(c23),PRE_ASSERT_(c24),\
PRE_ASSERT_(c25),PRE_ASSERT_(c26),PRE_ASSERT_(c27),PRE_ASSERT_(c28),\
PRE_ASSERT_(c29),PRE_ASSERT_(c30), &unused_assert_!=0)



#endif /* ! __cplusplus */



#define POST(c) EVALUATE_INVARIANT,POST0(c)
#define POST2(c1,c2) EVALUATE_INVARIANT,POST02(c1,c2)
#define POST3(c1,c2,c3) EVALUATE_INVARIANT,POST03(c1,c2,c3)
#define POST4(c1,c2,c3,c4) EVALUATE_INVARIANT,POST04(c1,c2,c3,c4)
#define POST5(c1,c2,c3,c4,c5) EVALUATE_INVARIANT,POST05(c1,c2,c3,c4,c5)
#define POST6(c1,c2,c3,c4,c5,c6) EVALUATE_INVARIANT,POST06(c1,c2,c3,c4,c5,c6)

#define POST7(c1,c2,c3,c4,c5,c6,c7) \
EVALUATE_INVARIANT,POST07(c1,c2,c3,c4,c5,c6,c7)

#define POST8(c1,c2,c3,c4,c5,c6,c7,c8) \
EVALUATE_INVARIANT,POST08(c1,c2,c3,c4,c5,c6,c7,c8)

#define POST9(c1,c2,c3,c4,c5,c6,c7,c8,c9) \
EVALUATE_INVARIANT,POST09(c1,c2,c3,c4,c5,c6,c7,c8,c9)

#define POST10(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10) \
EVALUATE_INVARIANT,POST010(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)

#define POST11(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11) \
EVALUATE_INVARIANT,POST011(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11)

#define POST12(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12) \
EVALUATE_INVARIANT,POST012(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12)

#define POST13(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13) \
EVALUATE_INVARIANT,POST013(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13)

#define POST14(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14) \
EVALUATE_INVARIANT,POST014(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14)

#define POST15(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15) \
EVALUATE_INVARIANT,POST015(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15)

#define POST16(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16) \
EVALUATE_INVARIANT,\
POST016(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16)

#define POST17(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17) \
EVALUATE_INVARIANT,\
POST017(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17)

#define POST18(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18) \
EVALUATE_INVARIANT,\
POST018(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18)

#define POST19(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19) \
EVALUATE_INVARIANT,\
POST019(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19)

#define POST20(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20) \
EVALUATE_INVARIANT,\
POST020(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20)

#define POST21(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21) \
EVALUATE_INVARIANT,\
POST021(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21)

#define POST22(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22) \
EVALUATE_INVARIANT,\
POST022(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22)

#define POST23(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23) \
EVALUATE_INVARIANT,\
POST023(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23)

#define POST24(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24) \
EVALUATE_INVARIANT,\
POST024(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24)

#define POST25(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25) \
EVALUATE_INVARIANT,\
POST025(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25)

#define POST26(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26) \
EVALUATE_INVARIANT,\
POST026(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26)

#define POST27(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27) \
EVALUATE_INVARIANT,\
POST027(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27)

#define POST28(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28) \
EVALUATE_INVARIANT,\
POST028(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28)

#define POST29(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29) \
EVALUATE_INVARIANT,\
POST029(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29)

#define POST30(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30) \
EVALUATE_INVARIANT,\
POST030(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,\
c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30)



#define CHECKING(s) s

#define OLD(variable) variable##_old_
#define REMEMBER(type,variable) type OLD(variable) = variable
#define REMEMBER_F(type,function_name) type OLD(function_name)=function_name()


#endif /* NO_ASSERTIONS */


/* Handle optional debugging statements separately. */

#ifdef DEBUGGING
#define DEBUG(s) s
#else
#define DEBUG(s)
#endif

#if defined( DEBUGGING ) || defined( DEBUGGING2 )
#define DEBUG2(s) s
#else
#define DEBUG2(unused_)
#endif



/******************************************************************************
 * Other macros useful in PRE() and POST() expressions:
 * E.g., PRE3( str, str[0], IN_RANGE( i, 0, 8 ) )
 *
 * Implies: p implies c: if p is true then c must be true.
 * Example:
 * void f( size_t count, int items[] ) { PRE( IMPLIES( count, items ) ) }
 */

#define IMPLIES(p,c) (!(p)||(c))
#define IMPLIES_ELSE(p,c1,c2) (((p)&&(c1))||((!(p))&&(c2)))

/* Boolean functions: */

#define NOT(a) (!(a))

/* multi-argument NOT() forms omitted due to ambiguous interpretation. */


#define AND2(a,b) ((a)&&(b))
#define AND3(a,b,c) ((a)&&(b)&&(c))
#define AND4(a,b,c,d) ((a)&&(b)&&(c)&&(d))
#define AND5(a,b,c,d,e) ((a)&&(b)&&(c)&&(d)&&(e))
#define AND6(a,b,c,d,e,f) ((a)&&(b)&&(c)&&(d)&&(e)&&(f))
#define AND7(a,b,c,d,e,f,g) ((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g))
#define AND8(a,b,c,d,e,f,g,h) ((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h))

#define AND9(a,b,c,d,e,f,g,h,i) \
((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i))

#define AND10(a,b,c,d,e,f,g,h,i,j) \
((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i)&&(j))

#define AND11(a,b,c,d,e,f,g,h,i,j,k) \
((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i)&&(j)&&(k))

#define AND12(a,b,c,d,e,f,g,h,i,j,k,l) \
((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i)&&(j)&&(k)&&(l))

#define AND13(a,b,c,d,e,f,g,h,i,j,k,l,m) \
((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i)&&(j)&&(k)&&(l)&&(m))

#define AND14(a,b,c,d,e,f,g,h,i,j,k,l,m,n) \
((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i)&&(j)&&(k)&&(l)&&(m)&&(n))

#define AND15(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) \
((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i)&&(j)&&(k)&&(l)&&(m)&&(n)&&(o))

#define AND16(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) \
((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i)&&(j)&&(k)&&(l)&&(m)&&(n)&&(o)&&(p))

#define AND17(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q) \
((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i)&&(j)&&(k)&&(l)&&(m)&&(n)&&(o)&&(p)\
&&(q))

#define AND18(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r) \
((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i)&&(j)&&(k)&&(l)&&(m)&&(n)&&(o)&&(p)\
&&(q)&&(r))

#define AND19(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s) \
((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i)&&(j)&&(k)&&(l)&&(m)&&(n)&&(o)&&(p)\
&&(q)&&(r)&&(s))

#define AND20(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t) \
((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i)&&(j)&&(k)&&(l)&&(m)&&(n)&&(o)&&(p)\
&&(q)&&(r)&&(s)&&(t))

#define AND21(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u) \
((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i)&&(j)&&(k)&&(l)&&(m)&&(n)&&(o)&&(p)\
&&(q)&&(r)&&(s)&&(t)&&(u))

#define AND22(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v) \
((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i)&&(j)&&(k)&&(l)&&(m)&&(n)&&(o)&&(p)\
&&(q)&&(r)&&(s)&&(t)&&(u)&&(v))

#define AND23(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w) \
((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i)&&(j)&&(k)&&(l)&&(m)&&(n)&&(o)&&(p)\
&&(q)&&(r)&&(s)&&(t)&&(u)&&(v)&&(w))

#define AND24(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x) \
((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i)&&(j)&&(k)&&(l)&&(m)&&(n)&&(o)&&(p)\
&&(q)&&(r)&&(s)&&(t)&&(u)&&(v)&&(w)&&(x))

#define AND25(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y) \
((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i)&&(j)&&(k)&&(l)&&(m)&&(n)&&(o)&&(p)\
&&(q)&&(r)&&(s)&&(t)&&(u)&&(v)&&(w)&&(x)&&(y))

#define AND26(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) \
((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i)&&(j)&&(k)&&(l)&&(m)&&(n)&&(o)&&(p)\
&&(q)&&(r)&&(s)&&(t)&&(u)&&(v)&&(w)&&(x)&&(y)&&(z))

#define AND27(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2) \
((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i)&&(j)&&(k)&&(l)&&(m)&&(n)&&(o)&&(p)\
&&(q)&&(r)&&(s)&&(t)&&(u)&&(v)&&(w)&&(x)&&(y)&&(z)&&(a2))

#define AND28(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,b2) \
((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i)&&(j)&&(k)&&(l)&&(m)&&(n)&&(o)&&(p)\
&&(q)&&(r)&&(s)&&(t)&&(u)&&(v)&&(w)&&(x)&&(y)&&(z)&&(a2)&&(b2))

#define AND29(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,b2,c2) \
((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i)&&(j)&&(k)&&(l)&&(m)&&(n)&&(o)&&(p)\
&&(q)&&(r)&&(s)&&(t)&&(u)&&(v)&&(w)&&(x)&&(y)&&(z)&&(a2)&&(b2)&&(c2))

#define AND30(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,b2,c2,d2) \
((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i)&&(j)&&(k)&&(l)&&(m)&&(n)&&(o)&&(p)\
&&(q)&&(r)&&(s)&&(t)&&(u)&&(v)&&(w)&&(x)&&(y)&&(z)&&(a2)&&(b2)&&(c2)&&(d2))



#define OR2(a,b) ((a)||(b))
#define OR3(a,b,c) ((a)||(b)||(c))
#define OR4(a,b,c,d) ((a)||(b)||(c)||(d))
#define OR5(a,b,c,d,e) ((a)||(b)||(c)||(d)||(e))
#define OR6(a,b,c,d,e,f) ((a)||(b)||(c)||(d)||(e)||(f))
#define OR7(a,b,c,d,e,f,g) ((a)||(b)||(c)||(d)||(e)||(f)||(g))
#define OR8(a,b,c,d,e,f,g,h) ((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h))

#define OR9(a,b,c,d,e,f,g,h,i) \
((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h)||(i))

#define OR10(a,b,c,d,e,f,g,h,i,j) \
((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h)||(i)||(j))

#define OR11(a,b,c,d,e,f,g,h,i,j,k) \
((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h)||(i)||(j)||(k))

#define OR12(a,b,c,d,e,f,g,h,i,j,k,l) \
((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h)||(i)||(j)||(k)||(l))

#define OR13(a,b,c,d,e,f,g,h,i,j,k,l,m) \
((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h)||(i)||(j)||(k)||(l)||(m))

#define OR14(a,b,c,d,e,f,g,h,i,j,k,l,m,n) \
((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h)||(i)||(j)||(k)||(l)||(m)||(n))

#define OR15(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) \
((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h)||(i)||(j)||(k)||(l)||(m)||(n)||(o))

#define OR16(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) \
((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h)||(i)||(j)||(k)||(l)||(m)||(n)||(o)||(p))

#define OR17(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q) \
((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h)||(i)||(j)||(k)||(l)||(m)||(n)||(o)||(p)\
||(q))

#define OR18(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r) \
((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h)||(i)||(j)||(k)||(l)||(m)||(n)||(o)||(p)\
||(q)||(r))

#define OR19(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s) \
((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h)||(i)||(j)||(k)||(l)||(m)||(n)||(o)||(p)\
||(q)||(r)||(s))

#define OR20(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t) \
((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h)||(i)||(j)||(k)||(l)||(m)||(n)||(o)||(p)\
||(q)||(r)||(s)||(t))

#define OR21(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u) \
((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h)||(i)||(j)||(k)||(l)||(m)||(n)||(o)||(p)\
||(q)||(r)||(s)||(t)||(u))

#define OR22(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v) \
((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h)||(i)||(j)||(k)||(l)||(m)||(n)||(o)||(p)\
||(q)||(r)||(s)||(t)||(u)||(v))

#define OR23(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w) \
((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h)||(i)||(j)||(k)||(l)||(m)||(n)||(o)||(p)\
||(q)||(r)||(s)||(t)||(u)||(v)||(w))

#define OR24(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x) \
((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h)||(i)||(j)||(k)||(l)||(m)||(n)||(o)||(p)\
||(q)||(r)||(s)||(t)||(u)||(v)||(w)||(x))

#define OR25(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y) \
((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h)||(i)||(j)||(k)||(l)||(m)||(n)||(o)||(p)\
||(q)||(r)||(s)||(t)||(u)||(v)||(w)||(x)||(y))

#define OR26(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) \
((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h)||(i)||(j)||(k)||(l)||(m)||(n)||(o)||(p)\
||(q)||(r)||(s)||(t)||(u)||(v)||(w)||(x)||(y)||(z))

#define OR27(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2) \
((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h)||(i)||(j)||(k)||(l)||(m)||(n)||(o)||(p)\
||(q)||(r)||(s)||(t)||(u)||(v)||(w)||(x)||(y)||(z)||(a2))

#define OR28(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,b2) \
((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h)||(i)||(j)||(k)||(l)||(m)||(n)||(o)||(p)\
||(q)||(r)||(s)||(t)||(u)||(v)||(w)||(x)||(y)||(z)||(a2)||(b2))

#define OR29(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,b2,c2) \
((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h)||(i)||(j)||(k)||(l)||(m)||(n)||(o)||(p)\
||(q)||(r)||(s)||(t)||(u)||(v)||(w)||(x)||(y)||(z)||(a2)||(b2)||(c2))

#define OR30(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,b2,c2,d2) \
((a)||(b)||(c)||(d)||(e)||(f)||(g)||(h)||(i)||(j)||(k)||(l)||(m)||(n)||(o)||(p)\
||(q)||(r)||(s)||(t)||(u)||(v)||(w)||(x)||(y)||(z)||(a2)||(b2)||(c2)||(d2))






/* Also see XOR() below. */


/* Numeric tests: */

/* Range membership: */

#define IN_RANGE(x,low,high) ((low)<=(x)&&(x)<=(high))
#define SIGN(x) ((x)<0?-1:1)

/* For C int 'flag' usage: */

#define IS_BOOL(a) ((a)==0||(a)==1)
#define IS_BOOL2(a,b) AND2(IS_BOOL(a),IS_BOOL(b))
#define IS_BOOL3(a,b,c) AND3(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c))
#define IS_BOOL4(a,b,c,d) AND4(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d))

#define IS_BOOL5(a,b,c,d,e) \
AND5(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e))

#define IS_BOOL6(a,b,c,d,e,f) \
AND6(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f))

#define IS_BOOL7(a,b,c,d,e,f,g) \
AND7(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g))

#define IS_BOOL8(a,b,c,d,e,f,g,h) \
AND8(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h))

#define IS_BOOL9(a,b,c,d,e,f,g,h,i) \
AND9(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h),IS_BOOL(i))

#define IS_BOOL10(a,b,c,d,e,f,g,h,i,j) \
AND10(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h),IS_BOOL(i),IS_BOOL(j))

#define IS_BOOL11(a,b,c,d,e,f,g,h,i,j,k) \
AND11(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h),IS_BOOL(i),IS_BOOL(j),IS_BOOL(k))

#define IS_BOOL12(a,b,c,d,e,f,g,h,i,j,k,l) \
AND12(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h),IS_BOOL(i),IS_BOOL(j),IS_BOOL(k),IS_BOOL(l))

#define IS_BOOL13(a,b,c,d,e,f,g,h,i,j,k,l,m) \
AND13(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h),IS_BOOL(i),IS_BOOL(j),IS_BOOL(k),IS_BOOL(l),IS_BOOL(m))

#define IS_BOOL14(a,b,c,d,e,f,g,h,i,j,k,l,m,n) \
AND14(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h),IS_BOOL(i),IS_BOOL(j),IS_BOOL(k),IS_BOOL(l),IS_BOOL(m),\
IS_BOOL(n))

#define IS_BOOL15(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) \
AND15(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h),IS_BOOL(i),IS_BOOL(j),IS_BOOL(k),IS_BOOL(l),IS_BOOL(m),\
IS_BOOL(n),IS_BOOL(o))

#define IS_BOOL16(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) \
AND16(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h),IS_BOOL(i),IS_BOOL(j),IS_BOOL(k),IS_BOOL(l),IS_BOOL(m),\
IS_BOOL(n),IS_BOOL(o),IS_BOOL(p))

#define IS_BOOL17(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q) \
AND17(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h),IS_BOOL(i),IS_BOOL(j),IS_BOOL(k),IS_BOOL(l),IS_BOOL(m),\
IS_BOOL(n),IS_BOOL(o),IS_BOOL(p),IS_BOOL(q))

#define IS_BOOL18(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r) \
AND18(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h),IS_BOOL(i),IS_BOOL(j),IS_BOOL(k),IS_BOOL(l),IS_BOOL(m),\
IS_BOOL(n),IS_BOOL(o),IS_BOOL(p),IS_BOOL(q),IS_BOOL(r))

#define IS_BOOL19(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s) \
AND19(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h),IS_BOOL(i),IS_BOOL(j),IS_BOOL(k),IS_BOOL(l),IS_BOOL(m),\
IS_BOOL(n),IS_BOOL(o),IS_BOOL(p),IS_BOOL(q),IS_BOOL(r),IS_BOOL(s))

#define IS_BOOL20(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t) \
AND20(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h),IS_BOOL(i),IS_BOOL(j),IS_BOOL(k),IS_BOOL(l),IS_BOOL(m),\
IS_BOOL(n),IS_BOOL(o),IS_BOOL(p),IS_BOOL(q),IS_BOOL(r),IS_BOOL(s),IS_BOOL(t))

#define IS_BOOL21(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u) \
AND21(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h),IS_BOOL(i),IS_BOOL(j),IS_BOOL(k),IS_BOOL(l),IS_BOOL(m),\
IS_BOOL(n),IS_BOOL(o),IS_BOOL(p),IS_BOOL(q),IS_BOOL(r),IS_BOOL(s),IS_BOOL(t),\
IS_BOOL(u))

#define IS_BOOL22(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v) \
AND22(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h),IS_BOOL(i),IS_BOOL(j),IS_BOOL(k),IS_BOOL(l),IS_BOOL(m),\
IS_BOOL(n),IS_BOOL(o),IS_BOOL(p),IS_BOOL(q),IS_BOOL(r),IS_BOOL(s),IS_BOOL(t),\
IS_BOOL(u),IS_BOOL(v))

#define IS_BOOL23(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w) \
AND23(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h),IS_BOOL(i),IS_BOOL(j),IS_BOOL(k),IS_BOOL(l),IS_BOOL(m),\
IS_BOOL(n),IS_BOOL(o),IS_BOOL(p),IS_BOOL(q),IS_BOOL(r),IS_BOOL(s),IS_BOOL(t),\
IS_BOOL(u),IS_BOOL(v),IS_BOOL(w))

#define IS_BOOL24(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x) \
AND24(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h),IS_BOOL(i),IS_BOOL(j),IS_BOOL(k),IS_BOOL(l),IS_BOOL(m),\
IS_BOOL(n),IS_BOOL(o),IS_BOOL(p),IS_BOOL(q),IS_BOOL(r),IS_BOOL(s),IS_BOOL(t),\
IS_BOOL(u),IS_BOOL(v),IS_BOOL(w),IS_BOOL(x))

#define IS_BOOL25(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y) \
AND25(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h),IS_BOOL(i),IS_BOOL(j),IS_BOOL(k),IS_BOOL(l),IS_BOOL(m),\
IS_BOOL(n),IS_BOOL(o),IS_BOOL(p),IS_BOOL(q),IS_BOOL(r),IS_BOOL(s),IS_BOOL(t),\
IS_BOOL(u),IS_BOOL(v),IS_BOOL(w),IS_BOOL(x),IS_BOOL(y))

#define IS_BOOL26(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) \
AND26(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h),IS_BOOL(i),IS_BOOL(j),IS_BOOL(k),IS_BOOL(l),IS_BOOL(m),\
IS_BOOL(n),IS_BOOL(o),IS_BOOL(p),IS_BOOL(q),IS_BOOL(r),IS_BOOL(s),IS_BOOL(t),\
IS_BOOL(u),IS_BOOL(v),IS_BOOL(w),IS_BOOL(x),IS_BOOL(y),IS_BOOL(z))

#define IS_BOOL27(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2) \
AND27(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h),IS_BOOL(i),IS_BOOL(j),IS_BOOL(k),IS_BOOL(l),IS_BOOL(m),\
IS_BOOL(n),IS_BOOL(o),IS_BOOL(p),IS_BOOL(q),IS_BOOL(r),IS_BOOL(s),IS_BOOL(t),\
IS_BOOL(u),IS_BOOL(v),IS_BOOL(w),IS_BOOL(x),IS_BOOL(y),IS_BOOL(z),IS_BOOL(a2))

#define IS_BOOL28(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,b2) \
AND28(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h),IS_BOOL(i),IS_BOOL(j),IS_BOOL(k),IS_BOOL(l),IS_BOOL(m),\
IS_BOOL(n),IS_BOOL(o),IS_BOOL(p),IS_BOOL(q),IS_BOOL(r),IS_BOOL(s),IS_BOOL(t),\
IS_BOOL(u),IS_BOOL(v),IS_BOOL(w),IS_BOOL(x),IS_BOOL(y),IS_BOOL(z),IS_BOOL(a2),\
IS_BOOL(b2))

#define IS_BOOL29(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,b2,\
c2)\
AND29(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h),IS_BOOL(i),IS_BOOL(j),IS_BOOL(k),IS_BOOL(l),IS_BOOL(m),\
IS_BOOL(n),IS_BOOL(o),IS_BOOL(p),IS_BOOL(q),IS_BOOL(r),IS_BOOL(s),IS_BOOL(t),\
IS_BOOL(u),IS_BOOL(v),IS_BOOL(w),IS_BOOL(x),IS_BOOL(y),IS_BOOL(z),IS_BOOL(a2),\
IS_BOOL(b2),IS_BOOL(c2))

#define IS_BOOL30(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,b2,\
c2,d2)\
AND30(IS_BOOL(a),IS_BOOL(b),IS_BOOL(c),IS_BOOL(d),IS_BOOL(e),IS_BOOL(f),\
IS_BOOL(g),IS_BOOL(h),IS_BOOL(i),IS_BOOL(j),IS_BOOL(k),IS_BOOL(l),IS_BOOL(m),\
IS_BOOL(n),IS_BOOL(o),IS_BOOL(p),IS_BOOL(q),IS_BOOL(r),IS_BOOL(s),IS_BOOL(t),\
IS_BOOL(u),IS_BOOL(v),IS_BOOL(w),IS_BOOL(x),IS_BOOL(y),IS_BOOL(z),IS_BOOL(a2),\
IS_BOOL(b2),IS_BOOL(c2),IS_BOOL(d2))



#define NON_ZERO2 AND2
#define NON_ZERO3 AND3
#define NON_ZERO4 AND4
#define NON_ZERO5 AND5
#define NON_ZERO6 AND6
#define NON_ZERO7 AND7
#define NON_ZERO8 AND8
#define NON_ZERO9 AND9
#define NON_ZERO10 AND10
#define NON_ZERO11 AND11
#define NON_ZERO12 AND12
#define NON_ZERO13 AND13
#define NON_ZERO14 AND14
#define NON_ZERO15 AND15
#define NON_ZERO16 AND16
#define NON_ZERO17 AND17
#define NON_ZERO18 AND18
#define NON_ZERO19 AND19
#define NON_ZERO20 AND20
#define NON_ZERO21 AND21
#define NON_ZERO22 AND22
#define NON_ZERO23 AND23
#define NON_ZERO24 AND24
#define NON_ZERO25 AND25
#define NON_ZERO26 AND26
#define NON_ZERO27 AND27
#define NON_ZERO28 AND28
#define NON_ZERO29 AND29
#define NON_ZERO30 AND30


#define IS_ZERO2(a,b) ((a)==0&&(b)==0)
#define IS_ZERO3(a,b,c) ((a)==0&&(b)==0&&(c)==0)
#define IS_ZERO4(a,b,c,d) ((a)==0&&(b)==0&&(c)==0&&(d)==0)
#define IS_ZERO5(a,b,c,d,e) ((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0)
#define IS_ZERO6(a,b,c,d,e,f) ((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0)

#define IS_ZERO7(a,b,c,d,e,f,g) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0)

#define IS_ZERO8(a,b,c,d,e,f,g,h) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0)

#define IS_ZERO9(a,b,c,d,e,f,g,h,i) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0&&(i)==0)

#define IS_ZERO10(a,b,c,d,e,f,g,h,i,j) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0&&(i)==0&&(j)==0)

#define IS_ZERO11(a,b,c,d,e,f,g,h,i,j,k) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0&&(i)==0&&\
(j)==0&&(k)==0)

#define IS_ZERO12(a,b,c,d,e,f,g,h,i,j,k,l) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0&&(i)==0&&\
(j)==0&&(k)==0&&(l)==0)

#define IS_ZERO13(a,b,c,d,e,f,g,h,i,j,k,l,m) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0&&(i)==0&&\
(j)==0&&(k)==0&&(l)==0&&(m)==0)

#define IS_ZERO14(a,b,c,d,e,f,g,h,i,j,k,l,m,n) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0&&(i)==0&&\
(j)==0&&(k)==0&&(l)==0&&(m)==0&&(n)==0)

#define IS_ZERO15(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0&&(i)==0&&\
(j)==0&&(k)==0&&(l)==0&&(m)==0&&(n)==0&&(o)==0)

#define IS_ZERO16(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0&&(i)==0&&\
(j)==0&&(k)==0&&(l)==0&&(m)==0&&(n)==0&&(o)==0&&(p)==0)

#define IS_ZERO17(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0&&(i)==0&&\
(j)==0&&(k)==0&&(l)==0&&(m)==0&&(n)==0&&(o)==0&&(p)==0&&(q)==0)

#define IS_ZERO18(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0&&(i)==0&&\
(j)==0&&(k)==0&&(l)==0&&(m)==0&&(n)==0&&(o)==0&&(p)==0&&(q)==0&&(r)==0)

#define IS_ZERO19(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0&&(i)==0&&\
(j)==0&&(k)==0&&(l)==0&&(m)==0&&(n)==0&&(o)==0&&(p)==0&&(q)==0&&(r)==0&&\
(s)==0)

#define IS_ZERO20(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0&&(i)==0&&\
(j)==0&&(k)==0&&(l)==0&&(m)==0&&(n)==0&&(o)==0&&(p)==0&&(q)==0&&(r)==0&&\
(s)==0&&(t)==0)

#define IS_ZERO21(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0&&(i)==0&&\
(j)==0&&(k)==0&&(l)==0&&(m)==0&&(n)==0&&(o)==0&&(p)==0&&(q)==0&&(r)==0&&\
(s)==0&&(t)==0&&(u)==0)

#define IS_ZERO22(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0&&(i)==0&&\
(j)==0&&(k)==0&&(l)==0&&(m)==0&&(n)==0&&(o)==0&&(p)==0&&(q)==0&&(r)==0&&\
(s)==0&&(t)==0&&(u)==0&&(v)==0)

#define IS_ZERO23(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0&&(i)==0&&\
(j)==0&&(k)==0&&(l)==0&&(m)==0&&(n)==0&&(o)==0&&(p)==0&&(q)==0&&(r)==0&&\
(s)==0&&(t)==0&&(u)==0&&(v)==0&&(w)==0)

#define IS_ZERO24(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0&&(i)==0&&\
(j)==0&&(k)==0&&(l)==0&&(m)==0&&(n)==0&&(o)==0&&(p)==0&&(q)==0&&(r)==0&&\
(s)==0&&(t)==0&&(u)==0&&(v)==0&&(w)==0&&(x)==0)

#define IS_ZERO25(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0&&(i)==0&&\
(j)==0&&(k)==0&&(l)==0&&(m)==0&&(n)==0&&(o)==0&&(p)==0&&(q)==0&&(r)==0&&\
(s)==0&&(t)==0&&(u)==0&&(v)==0&&(w)==0&&(x)==0&&(y)==0)

#define IS_ZERO26(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0&&(i)==0&&\
(j)==0&&(k)==0&&(l)==0&&(m)==0&&(n)==0&&(o)==0&&(p)==0&&(q)==0&&(r)==0&&\
(s)==0&&(t)==0&&(u)==0&&(v)==0&&(w)==0&&(x)==0&&(y)==0&&(z)==0)

#define IS_ZERO27(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0&&(i)==0&&\
(j)==0&&(k)==0&&(l)==0&&(m)==0&&(n)==0&&(o)==0&&(p)==0&&(q)==0&&(r)==0&&\
(s)==0&&(t)==0&&(u)==0&&(v)==0&&(w)==0&&(x)==0&&(y)==0&&(z)==0&&(a2)==0)

#define IS_ZERO28(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,b2) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0&&(i)==0&&\
(j)==0&&(k)==0&&(l)==0&&(m)==0&&(n)==0&&(o)==0&&(p)==0&&(q)==0&&(r)==0&&\
(s)==0&&(t)==0&&(u)==0&&(v)==0&&(w)==0&&(x)==0&&(y)==0&&(z)==0&&(a2)==0&&\
(b2)==0)

#define IS_ZERO29(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,b2,\
c2) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0&&(i)==0&&\
(j)==0&&(k)==0&&(l)==0&&(m)==0&&(n)==0&&(o)==0&&(p)==0&&(q)==0&&(r)==0&&\
(s)==0&&(t)==0&&(u)==0&&(v)==0&&(w)==0&&(x)==0&&(y)==0&&(z)==0&&(a2)==0&&\
(b2)==0&&(c2)==0)

#define IS_ZERO30(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,b2,\
c2,d2) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0&&(i)==0&&\
(j)==0&&(k)==0&&(l)==0&&(m)==0&&(n)==0&&(o)==0&&(p)==0&&(q)==0&&(r)==0&&\
(s)==0&&(t)==0&&(u)==0&&(v)==0&&(w)==0&&(x)==0&&(y)==0&&(z)==0&&(a2)==0&&\
(b2)==0&&(c2)==0&&(d2)==0)


#define GT_ZERO2(a,b) ((a)>0&&(b)>0)
#define GT_ZERO3(a,b,c) ((a)>0&&(b)>0&&(c)>0)
#define GT_ZERO4(a,b,c,d) ((a)>0&&(b)>0&&(c)>0&&(d)>0)
#define GT_ZERO5(a,b,c,d,e) ((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0)
#define GT_ZERO6(a,b,c,d,e,f) ((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0)

#define GT_ZERO7(a,b,c,d,e,f,g) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0)

#define GT_ZERO8(a,b,c,d,e,f,g,h) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0)

#define GT_ZERO9(a,b,c,d,e,f,g,h,i) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0&&(i)>0)

#define GT_ZERO10(a,b,c,d,e,f,g,h,i,j) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0&&(i)>0&&(j)>0)

#define GT_ZERO11(a,b,c,d,e,f,g,h,i,j,k) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0&&(i)>0&&(j)>0&&\
(k)>0)

#define GT_ZERO12(a,b,c,d,e,f,g,h,i,j,k,l) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0&&(i)>0&&(j)>0&&\
(k)>0&&(l)>0)

#define GT_ZERO13(a,b,c,d,e,f,g,h,i,j,k,l,m) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0&&(i)>0&&(j)>0&&\
(k)>0&&(l)>0&&(m)>0)

#define GT_ZERO14(a,b,c,d,e,f,g,h,i,j,k,l,m,n) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0&&(i)>0&&(j)>0&&\
(k)>0&&(l)>0&&(m)>0&&(n)>0)

#define GT_ZERO15(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0&&(i)>0&&(j)>0&&\
(k)>0&&(l)>0&&(m)>0&&(n)>0&&(o)>0)

#define GT_ZERO16(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0&&(i)>0&&(j)>0&&\
(k)>0&&(l)>0&&(m)>0&&(n)>0&&(o)>0&&(p)>0)

#define GT_ZERO17(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0&&(i)>0&&(j)>0&&\
(k)>0&&(l)>0&&(m)>0&&(n)>0&&(o)>0&&(p)>0&&(q)>0)

#define GT_ZERO18(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0&&(i)>0&&(j)>0&&\
(k)>0&&(l)>0&&(m)>0&&(n)>0&&(o)>0&&(p)>0&&(q)>0&&(r)>0)

#define GT_ZERO19(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0&&(i)>0&&(j)>0&&\
(k)>0&&(l)>0&&(m)>0&&(n)>0&&(o)>0&&(p)>0&&(q)>0&&(r)>0&&(s)>0)

#define GT_ZERO20(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0&&(i)>0&&(j)>0&&\
(k)>0&&(l)>0&&(m)>0&&(n)>0&&(o)>0&&(p)>0&&(q)>0&&(r)>0&&(s)>0&&(t)>0)

#define GT_ZERO21(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0&&(i)>0&&(j)>0&&\
(k)>0&&(l)>0&&(m)>0&&(n)>0&&(o)>0&&(p)>0&&(q)>0&&(r)>0&&(s)>0&&(t)>0&&(u)>0)

#define GT_ZERO22(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0&&(i)>0&&(j)>0&&\
(k)>0&&(l)>0&&(m)>0&&(n)>0&&(o)>0&&(p)>0&&(q)>0&&(r)>0&&(s)>0&&(t)>0&&(u)>0&&\
(v)>0)

#define GT_ZERO23(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0&&(i)>0&&(j)>0&&\
(k)>0&&(l)>0&&(m)>0&&(n)>0&&(o)>0&&(p)>0&&(q)>0&&(r)>0&&(s)>0&&(t)>0&&(u)>0&&\
(v)>0&&(w)>0)

#define GT_ZERO24(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0&&(i)>0&&(j)>0&&\
(k)>0&&(l)>0&&(m)>0&&(n)>0&&(o)>0&&(p)>0&&(q)>0&&(r)>0&&(s)>0&&(t)>0&&(u)>0&&\
(v)>0&&(w)>0&&(x)>0)

#define GT_ZERO25(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0&&(i)>0&&(j)>0&&\
(k)>0&&(l)>0&&(m)>0&&(n)>0&&(o)>0&&(p)>0&&(q)>0&&(r)>0&&(s)>0&&(t)>0&&(u)>0&&\
(v)>0&&(w)>0&&(x)>0&&(y)>0)

#define GT_ZERO26(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0&&(i)>0&&(j)>0&&\
(k)>0&&(l)>0&&(m)>0&&(n)>0&&(o)>0&&(p)>0&&(q)>0&&(r)>0&&(s)>0&&(t)>0&&(u)>0&&\
(v)>0&&(w)>0&&(x)>0&&(y)>0&&(z)>0)

#define GT_ZERO27(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0&&(i)>0&&(j)>0&&\
(k)>0&&(l)>0&&(m)>0&&(n)>0&&(o)>0&&(p)>0&&(q)>0&&(r)>0&&(s)>0&&(t)>0&&(u)>0&&\
(v)>0&&(w)>0&&(x)>0&&(y)>0&&(z)>0&&(a2)>0)

#define GT_ZERO28(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,b2) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0&&(i)>0&&(j)>0&&\
(k)>0&&(l)>0&&(m)>0&&(n)>0&&(o)>0&&(p)>0&&(q)>0&&(r)>0&&(s)>0&&(t)>0&&(u)>0&&\
(v)>0&&(w)>0&&(x)>0&&(y)>0&&(z)>0&&(a2)>0&&(b2)>0)

#define GT_ZERO29(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,b2,\
c2) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0&&(i)>0&&(j)>0&&\
(k)>0&&(l)>0&&(m)>0&&(n)>0&&(o)>0&&(p)>0&&(q)>0&&(r)>0&&(s)>0&&(t)>0&&(u)>0&&\
(v)>0&&(w)>0&&(x)>0&&(y)>0&&(z)>0&&(a2)>0&&(b2)>0&&(c2)>0)

#define GT_ZERO30(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,b2,\
c2,d2) \
((a)>0&&(b)>0&&(c)>0&&(d)>0&&(e)>0&&(f)>0&&(g)>0&&(h)>0&&(i)>0&&(j)>0&&\
(k)>0&&(l)>0&&(m)>0&&(n)>0&&(o)>0&&(p)>0&&(q)>0&&(r)>0&&(s)>0&&(t)>0&&(u)>0&&\
(v)>0&&(w)>0&&(x)>0&&(y)>0&&(z)>0&&(a2)>0&&(b2)>0&&(c2)>0&&(d2)>0)




#define GE_ZERO2(a,b) ((a)>=0&&(b)>=0)
#define GE_ZERO3(a,b,c) ((a)>=0&&(b)>=0&&(c)>=0)
#define GE_ZERO4(a,b,c,d) ((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0)
#define GE_ZERO5(a,b,c,d,e) ((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0)
#define GE_ZERO6(a,b,c,d,e,f) ((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0)

#define GE_ZERO7(a,b,c,d,e,f,g) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0)

#define GE_ZERO8(a,b,c,d,e,f,g,h) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0)

#define GE_ZERO9(a,b,c,d,e,f,g,h,i) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0&&(i)>=0)

#define GE_ZERO10(a,b,c,d,e,f,g,h,i,j) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0&&(i)>=0&&(j)>=0)

#define GE_ZERO11(a,b,c,d,e,f,g,h,i,j,k) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0&&(i)>=0&&\
(j)>=0&&(k)>=0)

#define GE_ZERO12(a,b,c,d,e,f,g,h,i,j,k,l) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0&&(i)>=0&&\
(j)>=0&&(k)>=0&&(l)>=0)

#define GE_ZERO13(a,b,c,d,e,f,g,h,i,j,k,l,m) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0&&(i)>=0&&\
(j)>=0&&(k)>=0&&(l)>=0&&(m)>=0)

#define GE_ZERO14(a,b,c,d,e,f,g,h,i,j,k,l,m,n) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0&&(i)>=0&&\
(j)>=0&&(k)>=0&&(l)>=0&&(m)>=0&&(n)>=0)

#define GE_ZERO15(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0&&(i)>=0&&\
(j)>=0&&(k)>=0&&(l)>=0&&(m)>=0&&(n)>=0&&(o)>=0)

#define GE_ZERO16(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0&&(i)>=0&&\
(j)>=0&&(k)>=0&&(l)>=0&&(m)>=0&&(n)>=0&&(o)>=0&&(p)>=0)

#define GE_ZERO17(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0&&(i)>=0&&\
(j)>=0&&(k)>=0&&(l)>=0&&(m)>=0&&(n)>=0&&(o)>=0&&(p)>=0&&(q)>=0)

#define GE_ZERO18(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0&&(i)>=0&&\
(j)>=0&&(k)>=0&&(l)>=0&&(m)>=0&&(n)>=0&&(o)>=0&&(p)>=0&&(q)>=0&&(r)>=0)

#define GE_ZERO19(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0&&(i)>=0&&\
(j)>=0&&(k)>=0&&(l)>=0&&(m)>=0&&(n)>=0&&(o)>=0&&(p)>=0&&(q)>=0&&(r)>=0&&(s)>=0)

#define GE_ZERO20(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0&&(i)>=0&&\
(j)>=0&&(k)>=0&&(l)>=0&&(m)>=0&&(n)>=0&&(o)>=0&&(p)>=0&&(q)>=0&&(r)>=0&&\
(s)>=0&&(t)>=0)

#define GE_ZERO21(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0&&(i)>=0&&\
(j)>=0&&(k)>=0&&(l)>=0&&(m)>=0&&(n)>=0&&(o)>=0&&(p)>=0&&(q)>=0&&(r)>=0&&\
(s)>=0&&(t)>=0&&(u)>=0)

#define GE_ZERO22(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0&&(i)>=0&&\
(j)>=0&&(k)>=0&&(l)>=0&&(m)>=0&&(n)>=0&&(o)>=0&&(p)>=0&&(q)>=0&&(r)>=0&&\
(s)>=0&&(t)>=0&&(u)>=0&&(v)>=0)

#define GE_ZERO23(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0&&(i)>=0&&\
(j)>=0&&(k)>=0&&(l)>=0&&(m)>=0&&(n)>=0&&(o)>=0&&(p)>=0&&(q)>=0&&(r)>=0&&\
(s)>=0&&(t)>=0&&(u)>=0&&(v)>=0&&(w)>=0)

#define GE_ZERO24(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0&&(i)>=0&&\
(j)>=0&&(k)>=0&&(l)>=0&&(m)>=0&&(n)>=0&&(o)>=0&&(p)>=0&&(q)>=0&&(r)>=0&&\
(s)>=0&&(t)>=0&&(u)>=0&&(v)>=0&&(w)>=0&&(x)>=0)

#define GE_ZERO25(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0&&(i)>=0&&\
(j)>=0&&(k)>=0&&(l)>=0&&(m)>=0&&(n)>=0&&(o)>=0&&(p)>=0&&(q)>=0&&(r)>=0&&\
(s)>=0&&(t)>=0&&(u)>=0&&(v)>=0&&(w)>=0&&(x)>=0&&(y)>=0)

#define GE_ZERO26(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0&&(i)>=0&&\
(j)>=0&&(k)>=0&&(l)>=0&&(m)>=0&&(n)>=0&&(o)>=0&&(p)>=0&&(q)>=0&&(r)>=0&&\
(s)>=0&&(t)>=0&&(u)>=0&&(v)>=0&&(w)>=0&&(x)>=0&&(y)>=0&&(z)>=0)

#define GE_ZERO27(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0&&(i)>=0&&\
(j)>=0&&(k)>=0&&(l)>=0&&(m)>=0&&(n)>=0&&(o)>=0&&(p)>=0&&(q)>=0&&(r)>=0&&\
(s)>=0&&(t)>=0&&(u)>=0&&(v)>=0&&(w)>=0&&(x)>=0&&(y)>=0&&(z)>=0&&(a2)>=0)

#define GE_ZERO28(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,b2) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0&&(i)>=0&&\
(j)>=0&&(k)>=0&&(l)>=0&&(m)>=0&&(n)>=0&&(o)>=0&&(p)>=0&&(q)>=0&&(r)>=0&&\
(s)>=0&&(t)>=0&&(u)>=0&&(v)>=0&&(w)>=0&&(x)>=0&&(y)>=0&&(z)>=0&&(a2)>=0&&\
(b2)>=0)

#define GE_ZERO29(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,b2,\
c2) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0&&(i)>=0&&\
(j)>=0&&(k)>=0&&(l)>=0&&(m)>=0&&(n)>=0&&(o)>=0&&(p)>=0&&(q)>=0&&(r)>=0&&\
(s)>=0&&(t)>=0&&(u)>=0&&(v)>=0&&(w)>=0&&(x)>=0&&(y)>=0&&(z)>=0&&(a2)>=0&&\
(b2)>=0&&(c2)>=0)

#define GE_ZERO30(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,b2,\
c2,d2) \
((a)>=0&&(b)>=0&&(c)>=0&&(d)>=0&&(e)>=0&&(f)>=0&&(g)>=0&&(h)>=0&&(i)>=0&&\
(j)>=0&&(k)>=0&&(l)>=0&&(m)>=0&&(n)>=0&&(o)>=0&&(p)>=0&&(q)>=0&&(r)>=0&&\
(s)>=0&&(t)>=0&&(u)>=0&&(v)>=0&&(w)>=0&&(x)>=0&&(y)>=0&&(z)>=0&&(a2)>=0&&\
(b2)>=0&&(c2)>=0&&(d2)>=0)


/* Set membership. */

#define IN3(x,a,b ) ((x)==(a)||(x)==(b))
#define IN4(x,a,b,c ) ((x)==(a)||(x)==(b)||(x)==(c))
#define IN5(x,a,b,c,d ) ((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d))

#define IN6(x,a,b,c,d,e) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e))

#define IN7(x,a,b,c,d,e,f) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f))

#define IN8(x,a,b,c,d,e,f,g) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g))

#define IN9(x,a,b,c,d,e,f,g,h) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||(x)==(h))

#define IN10(x,a,b,c,d,e,f,g,h,i) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||(x)==(h)\
||(x)==(i))

#define IN11(x,a,b,c,d,e,f,g,h,i,j) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||(x)==(h)\
||(x)==(i)||(x)==(j))

#define IN12(x,a,b,c,d,e,f,g,h,i,j,k) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||(x)==(h)\
||(x)==(i)||(x)==(j)||(x)==(k))

#define IN13(x,a,b,c,d,e,f,g,h,i,j,k,l) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||(x)==(h)\
||(x)==(i)||(x)==(j)||(x)==(k)||(x)==(l))

#define IN14(x,a,b,c,d,e,f,g,h,i,j,k,l,m) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||(x)==(h)\
||(x)==(i)||(x)==(j)||(x)==(k)||(x)==(l)||(x)==(m))

#define IN15(x,a,b,c,d,e,f,g,h,i,j,k,l,m,n) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||(x)==(h)\
||(x)==(i)||(x)==(j)||(x)==(k)||(x)==(l)||(x)==(m)||(x)==(n))

#define IN16(x,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||(x)==(h)\
||(x)==(i)||(x)==(j)||(x)==(k)||(x)==(l)||(x)==(m)||(x)==(n)||(x)==(o))

#define IN17(x,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||(x)==(h)\
||(x)==(i)||(x)==(j)||(x)==(k)||(x)==(l)||(x)==(m)||(x)==(n)||(x)==(o)\
||(x)==(p))

#define IN18(x,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||(x)==(h)\
||(x)==(i)||(x)==(j)||(x)==(k)||(x)==(l)||(x)==(m)||(x)==(n)||(x)==(o)\
||(x)==(p)||(x)==(q))

#define IN19(x,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||(x)==(h)\
||(x)==(i)||(x)==(j)||(x)==(k)||(x)==(l)||(x)==(m)||(x)==(n)||(x)==(o)\
||(x)==(p)||(x)==(q)||(x)==(r))

#define IN20(x,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||(x)==(h)\
||(x)==(i)||(x)==(j)||(x)==(k)||(x)==(l)||(x)==(m)||(x)==(n)||(x)==(o)\
||(x)==(p)||(x)==(q)||(x)==(r)||(x)==(s))

#define IN21(x,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||(x)==(h)\
||(x)==(i)||(x)==(j)||(x)==(k)||(x)==(l)||(x)==(m)||(x)==(n)||(x)==(o)\
||(x)==(p)||(x)==(q)||(x)==(r)||(x)==(s)||(x)==(t))

#define IN22(x,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||(x)==(h)\
||(x)==(i)||(x)==(j)||(x)==(k)||(x)==(l)||(x)==(m)||(x)==(n)||(x)==(o)\
||(x)==(p)||(x)==(q)||(x)==(r)||(x)==(s)||(x)==(t)||(x)==(u))

#define IN23(x,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||(x)==(h)\
||(x)==(i)||(x)==(j)||(x)==(k)||(x)==(l)||(x)==(m)||(x)==(n)||(x)==(o)\
||(x)==(p)||(x)==(q)||(x)==(r)||(x)==(s)||(x)==(t)||(x)==(u)||(x)==(v))

#define IN24(x,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||(x)==(h)\
||(x)==(i)||(x)==(j)||(x)==(k)||(x)==(l)||(x)==(m)||(x)==(n)||(x)==(o)\
||(x)==(p)||(x)==(q)||(x)==(r)||(x)==(s)||(x)==(t)||(x)==(u)||(x)==(v)\
||(x)==(w))

#define IN25(x,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,a2) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||(x)==(h)\
||(x)==(i)||(x)==(j)||(x)==(k)||(x)==(l)||(x)==(m)||(x)==(n)||(x)==(o)\
||(x)==(p)||(x)==(q)||(x)==(r)||(x)==(s)||(x)==(t)||(x)==(u)||(x)==(v)\
||(x)==(w)||(x)==(a2))

#define IN26(x,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,a2,b2) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||(x)==(h)\
||(x)==(i)||(x)==(j)||(x)==(k)||(x)==(l)||(x)==(m)||(x)==(n)||(x)==(o)\
||(x)==(p)||(x)==(q)||(x)==(r)||(x)==(s)||(x)==(t)||(x)==(u)||(x)==(v)\
||(x)==(w)||(x)==(a2)||(x)==(b2))

#define IN27(x,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,a2,b2,c2) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||(x)==(h)\
||(x)==(i)||(x)==(j)||(x)==(k)||(x)==(l)||(x)==(m)||(x)==(n)||(x)==(o)\
||(x)==(p)||(x)==(q)||(x)==(r)||(x)==(s)||(x)==(t)||(x)==(u)||(x)==(v)\
||(x)==(w)||(x)==(a2)||(x)==(b2)||(x)==(c2))

#define IN28(x,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,a2,b2,c2,d2) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||(x)==(h)\
||(x)==(i)||(x)==(j)||(x)==(k)||(x)==(l)||(x)==(m)||(x)==(n)||(x)==(o)\
||(x)==(p)||(x)==(q)||(x)==(r)||(x)==(s)||(x)==(t)||(x)==(u)||(x)==(v)\
||(x)==(w)||(x)==(a2)||(x)==(b2)||(x)==(c2)||(x)==(d2))

#define IN29(x,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,a2,b2,c2,d2,e2) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||(x)==(h)\
||(x)==(i)||(x)==(j)||(x)==(k)||(x)==(l)||(x)==(m)||(x)==(n)||(x)==(o)\
||(x)==(p)||(x)==(q)||(x)==(r)||(x)==(s)||(x)==(t)||(x)==(u)||(x)==(v)\
||(x)==(w)||(x)==(a2)||(x)==(b2)||(x)==(c2)||(x)==(d2)||(x)==(e2))

#define IN30(x,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,a2,b2,c2,d2,e2,f2)\
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||(x)==(h)\
||(x)==(i)||(x)==(j)||(x)==(k)||(x)==(l)||(x)==(m)||(x)==(n)||(x)==(o)\
||(x)==(p)||(x)==(q)||(x)==(r)||(x)==(s)||(x)==(t)||(x)==(u)||(x)==(v)\
||(x)==(w)||(x)==(a2)||(x)==(b2)||(x)==(c2)||(x)==(d2)||(x)==(e2)||(x)==(f2))


/* NON_ZERO_COUNT*() evaluates to the number of arguments that are non-zero: */

#define NON_ZERO_COUNT2(a,b) (((a)!=0)+((b)!=0))
#define NON_ZERO_COUNT3(a,b,c) (((a)!=0)+((b)!=0)+((c)!=0))
#define NON_ZERO_COUNT4(a,b,c,d) (((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0))

#define NON_ZERO_COUNT5(a,b,c,d,e) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0))

#define NON_ZERO_COUNT6(a,b,c,d,e,f) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0))

#define NON_ZERO_COUNT7(a,b,c,d,e,f,g) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0))

#define NON_ZERO_COUNT8(a,b,c,d,e,f,g,h) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0))

#define NON_ZERO_COUNT9(a,b,c,d,e,f,g,h,i) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0)\
+((i)!=0))

#define NON_ZERO_COUNT10(a,b,c,d,e,f,g,h,i,j) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0)\
+((i)!=0)+((j)!=0))

#define NON_ZERO_COUNT11(a,b,c,d,e,f,g,h,i,j,k) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0)\
+((i)!=0)+((j)!=0)+((k)!=0))

#define NON_ZERO_COUNT12(a,b,c,d,e,f,g,h,i,j,k,l) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0)\
+((i)!=0)+((j)!=0)+((k)!=0)+((l)!=0))

#define NON_ZERO_COUNT13(a,b,c,d,e,f,g,h,i,j,k,l,m) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0)\
+((i)!=0)+((j)!=0)+((k)!=0)+((l)!=0)+((m)!=0))

#define NON_ZERO_COUNT14(a,b,c,d,e,f,g,h,i,j,k,l,m,n) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0)\
+((i)!=0)+((j)!=0)+((k)!=0)+((l)!=0)+((m)!=0)+((n)!=0))

#define NON_ZERO_COUNT15(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0)\
+((i)!=0)+((j)!=0)+((k)!=0)+((l)!=0)+((m)!=0)+((n)!=0)+((o)!=0))

#define NON_ZERO_COUNT16(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0)\
+((i)!=0)+((j)!=0)+((k)!=0)+((l)!=0)+((m)!=0)+((n)!=0)+((o)!=0)+((p)!=0))

#define NON_ZERO_COUNT17(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0)\
+((i)!=0)+((j)!=0)+((k)!=0)+((l)!=0)+((m)!=0)+((n)!=0)+((o)!=0)+((p)!=0)\
+((q)!=0))

#define NON_ZERO_COUNT18(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0)\
+((i)!=0)+((j)!=0)+((k)!=0)+((l)!=0)+((m)!=0)+((n)!=0)+((o)!=0)+((p)!=0)\
+((q)!=0)+((r)!=0))

#define NON_ZERO_COUNT19(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0)\
+((i)!=0)+((j)!=0)+((k)!=0)+((l)!=0)+((m)!=0)+((n)!=0)+((o)!=0)+((p)!=0)\
+((q)!=0)+((r)!=0)+((s)!=0))

#define NON_ZERO_COUNT20(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0)\
+((i)!=0)+((j)!=0)+((k)!=0)+((l)!=0)+((m)!=0)+((n)!=0)+((o)!=0)+((p)!=0)\
+((q)!=0)+((r)!=0)+((s)!=0)+((t)!=0))

#define NON_ZERO_COUNT21(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0)\
+((i)!=0)+((j)!=0)+((k)!=0)+((l)!=0)+((m)!=0)+((n)!=0)+((o)!=0)+((p)!=0)\
+((q)!=0)+((r)!=0)+((s)!=0)+((t)!=0)+((u)!=0))

#define NON_ZERO_COUNT22(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0)\
+((i)!=0)+((j)!=0)+((k)!=0)+((l)!=0)+((m)!=0)+((n)!=0)+((o)!=0)+((p)!=0)\
+((q)!=0)+((r)!=0)+((s)!=0)+((t)!=0)+((u)!=0)+((v)!=0))

#define NON_ZERO_COUNT23(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0)\
+((i)!=0)+((j)!=0)+((k)!=0)+((l)!=0)+((m)!=0)+((n)!=0)+((o)!=0)+((p)!=0)\
+((q)!=0)+((r)!=0)+((s)!=0)+((t)!=0)+((u)!=0)+((v)!=0)+((w)!=0))

#define NON_ZERO_COUNT24(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0)\
+((i)!=0)+((j)!=0)+((k)!=0)+((l)!=0)+((m)!=0)+((n)!=0)+((o)!=0)+((p)!=0)\
+((q)!=0)+((r)!=0)+((s)!=0)+((t)!=0)+((u)!=0)+((v)!=0)+((w)!=0)+((x)!=0))

#define NON_ZERO_COUNT25(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0)\
+((i)!=0)+((j)!=0)+((k)!=0)+((l)!=0)+((m)!=0)+((n)!=0)+((o)!=0)+((p)!=0)\
+((q)!=0)+((r)!=0)+((s)!=0)+((t)!=0)+((u)!=0)+((v)!=0)+((w)!=0)+((x)!=0)\
+((y)!=0))

#define NON_ZERO_COUNT26(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0)\
+((i)!=0)+((j)!=0)+((k)!=0)+((l)!=0)+((m)!=0)+((n)!=0)+((o)!=0)+((p)!=0)\
+((q)!=0)+((r)!=0)+((s)!=0)+((t)!=0)+((u)!=0)+((v)!=0)+((w)!=0)+((x)!=0)\
+((y)!=0)+((z)!=0))

#define NON_ZERO_COUNT27(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
a2) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0)\
+((i)!=0)+((j)!=0)+((k)!=0)+((l)!=0)+((m)!=0)+((n)!=0)+((o)!=0)+((p)!=0)\
+((q)!=0)+((r)!=0)+((s)!=0)+((t)!=0)+((u)!=0)+((v)!=0)+((w)!=0)+((x)!=0)\
+((y)!=0)+((z)!=0)+((a2)!=0))

#define NON_ZERO_COUNT28(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
a2,b2) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0)\
+((i)!=0)+((j)!=0)+((k)!=0)+((l)!=0)+((m)!=0)+((n)!=0)+((o)!=0)+((p)!=0)\
+((q)!=0)+((r)!=0)+((s)!=0)+((t)!=0)+((u)!=0)+((v)!=0)+((w)!=0)+((x)!=0)\
+((y)!=0)+((z)!=0)+((a2)!=0)+((b2)!=0))

#define NON_ZERO_COUNT29(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
a2,b2,c2) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0)\
+((i)!=0)+((j)!=0)+((k)!=0)+((l)!=0)+((m)!=0)+((n)!=0)+((o)!=0)+((p)!=0)\
+((q)!=0)+((r)!=0)+((s)!=0)+((t)!=0)+((u)!=0)+((v)!=0)+((w)!=0)+((x)!=0)\
+((y)!=0)+((z)!=0)+((a2)!=0)+((b2)!=0)+((c2)!=0))

#define NON_ZERO_COUNT30(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,\
a2,b2,c2,d2) \
(((a)!=0)+((b)!=0)+((c)!=0)+((d)!=0)+((e)!=0)+((f)!=0)+((g)!=0)+((h)!=0)\
+((i)!=0)+((j)!=0)+((k)!=0)+((l)!=0)+((m)!=0)+((n)!=0)+((o)!=0)+((p)!=0)\
+((q)!=0)+((r)!=0)+((s)!=0)+((t)!=0)+((u)!=0)+((v)!=0)+((w)!=0)+((x)!=0)\
+((y)!=0)+((z)!=0)+((a2)!=0)+((b2)!=0)+((c2)!=0)+((d2)!=0))



/* Exclusive-or does not do short-circuited boolean evaluation: */

#define XOR2(a,b) (NON_ZERO_COUNT2(a,b)==1)
#define XOR3(a,b,c) (NON_ZERO_COUNT3(a,b,c)==1)
#define XOR4(a,b,c,d) (NON_ZERO_COUNT4(a,b,c,d)==1)
#define XOR5(a,b,c,d,e) (NON_ZERO_COUNT5(a,b,c,d,e)==1)
#define XOR6(a,b,c,d,e,f) (NON_ZERO_COUNT6(a,b,c,d,e,f)==1)
#define XOR7(a,b,c,d,e,f,g) (NON_ZERO_COUNT7(a,b,c,d,e,f,g)==1)
#define XOR8(a,b,c,d,e,f,g,h) (NON_ZERO_COUNT8(a,b,c,d,e,f,g,h)==1)
#define XOR9(a,b,c,d,e,f,g,h,i) (NON_ZERO_COUNT9(a,b,c,d,e,f,g,h,i)==1)
#define XOR10(a,b,c,d,e,f,g,h,i,j) (NON_ZERO_COUNT10(a,b,c,d,e,f,g,h,i,j)==1)

#define XOR11(a,b,c,d,e,f,g,h,i,j,k) \
(NON_ZERO_COUNT11(a,b,c,d,e,f,g,h,i,j,k)==1)

#define XOR12(a,b,c,d,e,f,g,h,i,j,k,l) \
(NON_ZERO_COUNT12(a,b,c,d,e,f,g,h,i,j,k,l)==1)

#define XOR13(a,b,c,d,e,f,g,h,i,j,k,l,m) \
(NON_ZERO_COUNT13(a,b,c,d,e,f,g,h,i,j,k,l,m)==1)

#define XOR14(a,b,c,d,e,f,g,h,i,j,k,l,m,n) \
(NON_ZERO_COUNT14(a,b,c,d,e,f,g,h,i,j,k,l,m,n)==1)

#define XOR15(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) \
(NON_ZERO_COUNT15(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o)==1)

#define XOR16(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) \
(NON_ZERO_COUNT16(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p)==1)

#define XOR17(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q) \
(NON_ZERO_COUNT17(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q)==1)

#define XOR18(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r) \
(NON_ZERO_COUNT18(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r)==1)

#define XOR19(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s) \
(NON_ZERO_COUNT19(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s)==1)

#define XOR20(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t) \
(NON_ZERO_COUNT20(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t)==1)

#define XOR21(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u) \
(NON_ZERO_COUNT21(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u)==1)

#define XOR22(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v) \
(NON_ZERO_COUNT22(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v)==1)

#define XOR23(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w) \
(NON_ZERO_COUNT23(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w)==1)

#define XOR24(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x) \
(NON_ZERO_COUNT24(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x)==1)

#define XOR25(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y) \
(NON_ZERO_COUNT25(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y)==1)

#define XOR26(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) \
(NON_ZERO_COUNT26(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z)==1)

#define XOR27(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2) \
(NON_ZERO_COUNT27(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2)==1)

#define XOR28(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,b2) \
(NON_ZERO_COUNT28(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,\
b2)==1)

#define XOR29(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,b2,c2) \
(NON_ZERO_COUNT29(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,\
b2,c2)==1)

#define XOR30(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,b2,c2,d2)\
(NON_ZERO_COUNT30(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a2,\
b2,c2,d2)==1)


#endif /* ASSERTIONS_H */



