import { NextRequest, NextResponse } from 'next/server'
import { getCachedPapers } from '@/lib/external-data-cache'

export async function GET(request: NextRequest) {
  try {
    const { searchParams } = new URL(request.url)
    const forceRefresh = searchParams.get('refresh') === 'true'

    const result = await getCachedPapers(forceRefresh)

    return NextResponse.json({
      papers: result.papers,
      total: result.total,
      fromCache: result.fromCache,
      timestamp: Date.now()
    })
  } catch (error) {
    console.error('Error in literature cache API:', error)

    // Return error response with empty data
    return NextResponse.json({
      papers: [],
      total: 0,
      fromCache: false,
      error: 'Failed to load literature data',
      timestamp: Date.now()
    }, { status: 500 })
  }
}