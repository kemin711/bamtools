#include "api/BamAlignment.h"
#include <gtest/gtest.h>

using namespace std;
using namespace BamTools;

class BamAlignmentTest : public testing::Test {
   public:
      BamAlignmentTest() {}
      virtual ~BamAlignmentTest() {}
      virtual void SetUp() { }
      virtual void TearDown() { }

      BamAlignment empty;
};

TEST_F(BamAlignmentTest, Subsequence) {
}

int main(int argc, char* argv[]) {
   testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
